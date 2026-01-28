#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import os                           # filesystem utilities
import pathlib                      # library for handle filesystem paths
import matplotlib.pyplot as plt     # library for visualization
import seaborn as sns               # library for statistical data visualization
import pandas as pd                 # library for data analysis and manipulation
import scanpy as sc                 # single-cell data processing
import scanpy.external as sce       # library for harmony integration
import muon as mu
import numpy as np
import anndata as ad
from matplotlib.backends.backend_pdf import PdfPages

warnings.filterwarnings("ignore")
# PARAMETERS

# set script version number
VERSION = "0.0.1"

# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function clusters single-cell data at different resolution levels.
    """

# --------------------------------------------------------------------------------------------------------------------
#                                          LIBRARY CONFIG
# --------------------------------------------------------------------------------------------------------------------

    sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

# Define command line arguments with argparse

    parser = argparse.ArgumentParser(description = "Clustering",
        epilog = "This function clusters cells at different resolutions and identifies marker genes for each cluster and resolution",
        )
    parser.add_argument('-ad','--input-h5mu-file', metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5MU_OUTPUT_FILE', type=pathlib.Path, default="matrix.clustered.h5mu",
                        help="name of the output h5ad file after clustering")
    parser.add_argument('-e', '--excel_out', metavar='RANKED_GENES_XLSX', default="ranked_genes.xlsx",
                        help="path and name of excel table with ranked marker genes for each cluster and resolution")
    parser.add_argument('-csv', '--csv_out', metavar='H5AD_OUTPUT_FILE', default="final_metadata.csv",
                        help="path and name of csv tabel with UMAP coordinates for each cell")
    parser.add_argument('-min_res', '--resolution_min',  dest='min_res', type=float, default=0.1, 
                        help="Minimum clustering resolution to be evaluated (default: 0.1). All resolution values between min_res and max_res will be evaluated.")
    parser.add_argument('-max_res', '--resolution_max',  dest='max_res', type=float, default=1.1, 
                        help="Maximum clustering resolution to be evaluated (default: 1.1). All resolution values between min_res and max_res will be evaluated.")
    parser.add_argument('-n', '--top_n',  dest='top_n', type=int, default=10, 
                        help="number of top marker genes to save for each cluster (default is 10)")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_file = args.input_h5mu_files
    output = args.out
    output_excel= args.excel_out
    output_csv= args.csv_out
    min_res = args.min_res
    max_res = args.max_res
    top_n = args.top_n

    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print("-File {}:".format(str(input_h5mu_file)))

# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5MU FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the MTX combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder ... ", end ='')
    mdata = mu.read_h5mu(input_h5mu_file)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    # Extracting gex, handling genes metadata and associated rownames. 
    # Double assignment is needed in order to keep both the annotations (symbols and ids)
    print("\n===== GEX MODALITY DATA =====")
    gex = mdata.mod['gex']
    gex.var['gene_id'] = gex.var.index.astype(str)
    gex.var_names = gex.var['gene_symbols'].astype(str)
    gex.var.index.name = None
    gex.var['gene_symbols'] = gex.var.index.astype(str)

# --------------------------------------------------------------------------------------------------------------------
#                                 CLUSTERING
# --------------------------------------------------------------------------------------------------------------------
    # Clustering was performed on the Harmony representation (integrated reduced matrix)
    print(gex.uns['neighbors']['params']['use_rep'])

    print("\n===== CLUSTERING =====")
    # Clusters cells based on transcriptional similarities
    print("\nComputing Leiden clustering at different resolutions")

    clustering_labels = []
    # Test all the resolution parameters between min_res and max_res
    resolutions = np.round(np.arange(min_res, max_res, 0.1), 2)
    
    for res in resolutions:
        clustering_labels.append("leiden_{}".format(res))
        if "leiden_{}".format(res) in gex.obs:
            print("leiden_{}".format(res) + " already exists... going on with next resolution.")
            continue
        sc.tl.leiden(gex, resolution=res, key_added="leiden_{}".format(res))

# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE MARKER GENES FOR EACH RESOLUTION AND FOR EACH CLUSTER
# --------------------------------------------------------------------------------------------------------------------

    heat_name = "top_" + str(top_n) + "_markers_heatmap.pdf"
    with pd.ExcelWriter(output_excel) as writer, PdfPages(os.path.join(args.results, heat_name)) as pdf_heatmap:
        for res in resolutions:
            
            #Compute marker genes for each cluster, expects logarithmized data
            print("\nComputing marker genes for each clusters at resolution {}".format(res))
            sc.tl.rank_genes_groups(gex, groupby="leiden_{}".format(res), method="wilcoxon", key_added="leiden_{}".format(res), pts=True)
            df = sc.get.rank_genes_groups_df(gex, group=None, pval_cutoff=0.05, log2fc_min=0.25, key="leiden_{}".format(res))
            #df['gene_symbols'] = df['names'].map(gex.var['gene_symbols'].to_dict())

            # Select the top 10 markers per cluster
            top = (df.sort_values(["group", "pvals_adj", "logfoldchanges"], ascending=[True, True, False]).groupby("group").head(top_n))
            top_genes = top["names"].tolist() # .unique().tolist()
            
            # Heatmap
            print("\nPlotting top marker genes heatmap for resolution {}".format(res))
            plt.figure(figsize=(45, 55))
            sc.pl.heatmap(gex, var_names=top_genes, groupby="leiden_{}".format(res), show=False, cmap='viridis')
            plt.title(f"Top {top_n} marker genes – resolution {res}")
            pdf_heatmap.savefig(bbox_inches="tight", dpi=300)
            plt.close()
            print("Done!")

            print("\nSaving marker genes for each cluster and resolution in excel file")
            df.to_excel(writer, sheet_name=f"Leiden_{res}", index=False)
            print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

    # Visualize Leiden clustering on UMAP plot
    print("\nVisualized Leiden clustering on UMAP plot")
    plot_name = "cluster_id_all.pdf"
    sc.pl.umap(gex, color=clustering_labels, legend_loc='on data', show=False)
    plt.savefig(os.path.join(args.results, plot_name), bbox_inches='tight', dpi=300)
    plt.close()

    plot_name_each = os.path.join(args.results, "cluster_id_each.pdf")
    with PdfPages(plot_name_each) as pdf:
        for res, cl in zip(resolutions, clustering_labels):
            print(f"\nVisualized Leiden clustering on UMAP plot at resolution {res}")
            plt.figure(figsize=(10, 8))
            sc.pl.umap(gex, color=cl, legend_loc='on data', show=False)
            pdf.savefig(bbox_inches='tight', dpi=300)
            plt.close()
        

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
    mdata.mod['gex'] = gex
    mdata.update()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    df = pd.DataFrame(gex.obsm["X_umap"], index=gex.obs_names).rename(columns={0: "X_UMAP", 1: "Y_UMAP"})
    df = df.join(gex.obs)
    df.index.name = 'cell_barcodes'
    print("Saving csv table with Harmony corrected UMAP coordinates for each cell and Leiden ID {}".format(output_csv))
    df.to_csv(output_csv)
    print("Done!")

    print("Saving h5ad data to file {}".format(output))
    mdata.write(output)
    print("Done!")


#####################################################################################################


if __name__ == '__main__':
    main()
