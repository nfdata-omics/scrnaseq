#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import os                           # filesystem utilities
import re                           # hanlding regex
import pathlib                      # library for handle filesystem paths
import matplotlib.pyplot as plt     # library for visualization
import pandas as pd                 # library for data analysis and manipulation
import scanpy as sc                 # single-cell data processing
import mudata as md
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
    This function perfors dimensionality reduction of dataset.
    """
# --------------------------------------------------------------------------------------------------------------------
#                                          LIBRARY CONFIG
# --------------------------------------------------------------------------------------------------------------------

    sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

# Define command line arguments with argparse

    parser = argparse.ArgumentParser(prog='DimRed', usage='%(prog)s [options]', description = "Feature selection and dimensionality reduction",
        epilog = "This function reduce the dimensionality of the dataset and only include the most informative genes.")
    parser.add_argument('-ad', '--input-h5mu-file', metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5MU_OUTPUT_FILE', type=pathlib.Path, default="matrix.hvg.h5mu",
                        help="name of the output h5ad file after dimensionality reduction")
    parser.add_argument('-csv', '--csv_out', metavar='CSV_TABLE',type=pathlib.Path, default="umap_coordinates.csv",
                        help="csv table with UMAP coordinates for each cell")
    parser.add_argument('-pcs', '--n_princ_comps', dest='n_pcs', type=int, default=30, help="number of principal components to select and use for UMAP algorithms")
    parser.add_argument('-nn', '--n_neighbors', dest='n_neighbors', type=int, default=20, help="Size of local neighborhood used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. Values should be in the range 2 to 100")
    parser.add_argument('-md', '--min_dist', dest='min_dist', type=float, default=0.1, help="minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together")
    parser.add_argument('-r', '--results', type=pathlib.Path, default=pathlib.Path('./'), help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_file = args.input_h5mu_files
    output = args.out
    output_csv= args.csv_out
    n_pcs = args.n_pcs
    n_neighbors = args.n_neighbors
    min_dist = args.min_dist

# print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print(f"-File {input_h5mu_file}:")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the MTX combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder ... ", end ='')
    mdata = md.read(input_h5mu_file)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== GEX MODALITY DATA =====")
    gex = mdata.mod['gex']
    gex.var["feature_types"].value_counts()

# --------------------------------------------------------------------------------------------------------------------
#                                 FEATURE SELECTION & DIMENSIONALITY REDUCTION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== FEATURE SELECTION =====")
    # Selecting genes that exhibit high variability consistently across different batches (i.e. samples)
    print("\nSelecting highly-variable genes")
    #sc.pp.highly_variable_genes(gex, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=False)
    try:
        top_10_percent_genes = max(int(0.1 * gex.shape[1]), 10)
        sc.pp.highly_variable_genes(gex, n_top_genes=top_10_percent_genes, batch_key="sample")
    except IndexError:
        gex.var['highly_variable'] = True

    num_hvg = gex.var['highly_variable'].sum()
    print("\nNumber of highly variable genes (top 10%) identified: ", num_hvg)

    print("\n===== DIMENSIONALITY REDUCTION =====")
    print("\nPerforming dimensionality reduction by running principal component analysis (PCA)")
    n_comps_to_compute = min(50, num_hvg - 1)
    sc.tl.pca(gex, use_highly_variable = True, n_comps=n_comps_to_compute)

    # Visualize PCA plot
    print("\nVisualized PCA plot")
    plt.figure(figsize=(35, 25))
    sc.pl.pca(gex, color='sample', show=False)
    plt.savefig(os.path.join(args.results,'pca_GEX.pdf'), bbox_inches='tight', dpi=300)
    plt.close()

    # Visualize Elbow plot
    print("\nVisualized Elbow plot for PCA components")
    plt.figure(figsize=(35, 25))
    sc.pl.pca_variance_ratio(gex, n_pcs=n_comps_to_compute, show=False)
    plt.savefig(os.path.join(args.results,'pca_elbow.pdf'), bbox_inches='tight', dpi=300)
    plt.close()
# --------------------------------------------------------------------------------------------------------------------
#                                 DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== NEAREST NEIGHBOR GRAPH CONSTRUCTION =====")
    print("\nConstruction of the nearest neighbor graph")
    sc.pp.neighbors(gex, n_neighbors=n_neighbors, n_pcs=n_pcs)

    print("\n===== DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION=====")
    print("\nPerforming dimensionality reduction by running uniform manifold approximation and projection (UMAP)")
    sc.tl.umap(gex, min_dist=min_dist, random_state=42)

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

    # Visualize UMAP plot
    print("\nVisualized UMAP plot")
    plt.figure(figsize=(35, 25))
    sc.pl.umap(gex, color ='sample', show=False)
    plt.savefig(os.path.join(args.results, 'umap_plot_GEX.pdf'), bbox_inches='tight', dpi=300)
    plt.close()

    # UMAP plot highlighting metadata features (if present)
    # Metadata features have been renamed as meta_* in the convert_mudata step
    pattern = re.compile(r"meta_.*")
    meta_cols = [col for col in gex.obs.columns if pattern.match(col)]
    if len(meta_cols) > 0:
        with PdfPages(os.path.join(args.results, "umap_plot_GEX_metadata.pdf")) as pdf:
            for col in meta_cols:
                plt.figure(figsize=(45, 35))
                sc.pl.umap(gex, color=col, show=False)
                pdf.savefig(bbox_inches="tight", dpi=300)
                plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
    mdata.mod['gex'] = gex
    mdata.update()

# --------------------------------------------------------------------------------------------------------------------
#                                 CITE MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    if 'pro' in mdata.mod:
        print("\n===== CITE MODALITY DATA =====")
        pro = mdata.mod['pro']
        pro.var["feature_types"].value_counts()

# --------------------------------------------------------------------------------------------------------------------
#                                 DIMENSIONALITY REDUCTION
# --------------------------------------------------------------------------------------------------------------------

        print("\n===== DIMENSIONALITY REDUCTION =====")
        print("\nPerforming dimensionality reduction by running principal component analysis (PCA)")
        sc.pp.pca(pro, svd_solver="arpack")
        sc.pl.pca_variance_ratio(pro, n_pcs=50)

# --------------------------------------------------------------------------------------------------------------------
#                                 DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION
# --------------------------------------------------------------------------------------------------------------------

        print("\n===== NEAREST NEIGHBOR GRAPH CONSTRUCTION =====")
        print("\nConstruction of the nearest neighbor graph")
        sc.pp.neighbors(pro, n_pcs=30)

        print("\n===== DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION=====")
        print("\nPerforming dimensionality reduction by running uniform manifold approximation and projection (UMAP)")
        sc.tl.umap(pro)


# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

        # Visualize UMAP plot
        print("\nVisualized UMAP plot")
        plt.figure(figsize=(14, 13))
        sc.pl.umap(pro, color ='sample', show=False, legend_loc='lower right', legend_fontsize=12)
        plt.savefig(os.path.join(args.results,'umap_plot_CITE.pdf'), bbox_inches='tight', dpi=300)
        plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE CITE DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
        print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
        mdata.mod['pro'] = pro
        mdata.update()
    else:
        print("CITE modality does not exist in mdata.mod.")

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    print(f"Saving h5mu data to file {output}")
    mdata.write(output)
    print("Done!")

    df = pd.DataFrame(gex.obsm["X_umap"], index=gex.obs_names).rename(columns={0: "X_UMAP", 1: "Y_UMAP"})
    df.index.name = 'cell_barcodes'
    print(f"Saving csv table with UMAP coordinates for each cell {output_csv}")
    df.to_csv(output_csv)
    print("Done!")


#####################################################################################################


if __name__ == '__main__':
    main()
