#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

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
    parser.add_argument('-ad','--input-h5mu-file',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5MU_OUTPUT_FILE', type=pathlib.Path, default="matrix.clustered.h5mu",
                        help="name of the output h5ad file after clustering")
    parser.add_argument('-e', '--excel_out', metavar='RANKED_GENES_XLSX', default="ranked_genes.xlsx",
                        help="path and name of excel table with ranked marker genes for each cluster and resolution")
    parser.add_argument('-csv', '--csv_out', metavar='H5AD_OUTPUT_FILE', default="final_metadata.csv",
                        help="path and name of csv tabel with UMAP coordinates for each cell")
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
    mdata= mu.read_h5mu(input_h5mu_file)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== GEX MODALITY DATA =====")
    gex = mdata.mod['gex']
# --------------------------------------------------------------------------------------------------------------------
#                                 CLUSTERING
# --------------------------------------------------------------------------------------------------------------------
    # Clustering was performed on the Harmony representation (integrated reduced matrix)
    print(gex.uns['neighbors']['params']['use_rep'])

    print("\n===== CLUSTERING =====")
    # Clusters cells based on transcriptional similarities
    print("\nComputing Leiden clustering at different resolutions")

    clustering_labels = []
    for res in np.round(np.arange(0.1, 1.0, 0.3),2):
        clustering_labels.append("leiden_{}".format(res))
        if "leiden_{}".format(res) in gex.obs:
            print("leiden_{}".format(res) + " already exists... going on with next resolution.")
            continue
        sc.tl.leiden(gex, resolution=res, key_added="leiden_{}".format(res))

# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE MARKER GENES FOR EACH RESOLUTION AND FOR EACH CLUSTER
# --------------------------------------------------------------------------------------------------------------------

    with pd.ExcelWriter(output_excel) as writer:

        for res in np.round(np.arange(0.1, 1.0, 0.1),2):
            print("\nComputing top 20 marker genes for each clusters at resolution {}".format(res))
            #Compute top 20 marker genes for each cluster, expects logarithmized data
            sc.tl.rank_genes_groups(gex, groupby="leiden_{}".format(res),method="wilcoxon",key_added="leiden_{}".format(res), n_genes=100,pts=True)
            df = sc.get.rank_genes_groups_df(gex, group=None,pval_cutoff=0.05, log2fc_min=0.25,key="leiden_{}".format(res))
            #df['gene_symbol'] = df['names'].map(gex.var['gene_symbols'].to_dict())

            print("\nSaving top 20 marker genes for each cluster and resolution in excel file")
            df.to_excel(writer, sheet_name=f"Leiden_{res}", index=False)
            print("Done!")


        print("\nComputing top 20 marker genes for each sample {}".format(res))
        #Compute top 20 marker genes for each sample, expects logarithmized data
        sc.tl.rank_genes_groups(gex, groupby="sample",method="wilcoxon",key_added="sample_marker", n_genes=100,pts=True)
        df = sc.get.rank_genes_groups_df(gex, group=None,pval_cutoff=0.05, log2fc_min=0.25,key="sample_marker")
        #df['gene_symbol'] = df['names'].map(gex.var['gene_symbols'].to_dict())

        print("\nSaving top 20 marker genes for each cluster and resolution in excel file")
        df.to_excel(writer, sheet_name="Sample_Markers",index=False)
        print("Done!")
# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

    # Visualize Leiden clustering on UMAP plot

    print("\nVisualized Leiden clustering on UMAP plot")
    sc.pl.umap(gex, color=clustering_labels ,legend_loc='on data',show=False)
    plt.savefig(os.path.join(args.results,'cluster_id.png'))
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
