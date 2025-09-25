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
import pandas as pd                 # library for data analysis and manipulation
import scanpy as sc                 # single-cell data processing
import mudata as md


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

    parser = argparse.ArgumentParser(prog='DimRed', usage='%(prog)s [options]',description = "Feature selection and dimensionality reduction",
        epilog = "This function reduce the dimensionality of the dataset and only include the most informative genes.")
    parser.add_argument('-ad','--input-h5mu-file',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5MU_OUTPUT_FILE', type=pathlib.Path, default="matrix.hvg.h5mu",
                        help="name of the output h5ad file after dimensionality reduction")
    parser.add_argument('-csv', '--csv_out', metavar='CSV_TABLE',type=pathlib.Path, default="umap_coordinates.csv",
                        help="csv tabel with UMAP coordinates for each cell")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_file = args.input_h5mu_files
    output = args.out
    output_csv=args.csv_out

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
    mdata= md.read(input_h5mu_file)
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
    # select highly=variable genes for each sample
    print("\nSelecting highly-variable genes are selected within each batch separately and merged")
    sc.pp.highly_variable_genes(gex, min_mean=0.0125, max_mean=3, min_disp=0.5,subset=False)

    print("\n===== DIMENSIONALITY REDUCTION =====")
    print("\nPerforming dimensionality reduction by running principal component analysis (PCA)")
    sc.tl.pca(gex,use_highly_variable = True, n_comps=50)

    # Visualize PCA plot
    print("\nVisualized PCA plot")
    plt.figure(figsize=(12, 10))
    sc.pl.pca(gex, color='sample',show=False)
    plt.savefig(os.path.join(args.results,'pca_GEX.png'))
    plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                                 DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== NEAREST NEIGHBOR GRAPH CONSTRUCTION =====")
    print("\nConstruction of the nearest neighbor graph")
    sc.pp.neighbors(gex,n_neighbors=30,n_pcs=20)

    print("\n===== DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION=====")
    print("\nPerforming dimensionality reduction by running uniform manifold approximation and projection (UMAP)")
    sc.tl.umap(gex,min_dist=0.1,random_state=42)

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

    # Visualize UMAP plot
    print("\nVisualized UMAP plot")
    plt.figure(figsize=(12, 10))
    sc.pl.umap(gex, color ='sample',show=False)
    plt.savefig(os.path.join(args.results,'umap_plot_GEX.png'))
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
        sc.pp.neighbors(pro, n_pcs=20)

        print("\n===== DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION=====")
        print("\nPerforming dimensionality reduction by running uniform manifold approximation and projection (UMAP)")
        sc.tl.umap(pro)


# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

        # Visualize UMAP plot
        print("\nVisualized UMAP plot")
        plt.figure(figsize=(12, 10))
        sc.pl.umap(pro, color ='sample',show=False)
        plt.savefig(os.path.join(args.results,'umap_plot_CITE.png'))
        plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE CITE DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
        print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
        mdata.mod['pro'] = pro
        mdata.update()
    else:
        print("CITE modality does not exist in mdata.mod.")

    '''
# --------------------------------------------------------------------------------------------------------------------
#                                 ATAC MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    if 'atac' in mdata.mod:
        print("\n===== ATAC MODALITY DATA =====")
        atac = mdata.mod['atac']

# --------------------------------------------------------------------------------------------------------------------
#                                 FEATURE SELECTION & DIMENSIONALITY REDUCTION
# --------------------------------------------------------------------------------------------------------------------

        print("\n===== FEATURE SELECTION =====")
        # select highly-variable peaks for each sample
        print("\nSelecting highly-variable paeks are selected within each batch separately and merged")
        sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.53, min_disp=5,batch_key = 'sample',subset=False)

        print("\n===== DIMENSIONALITY REDUCTION =====")
        print("\nPerforming dimensionality reduction by running principal component analysis (PCA)")
        sc.pp.scale(atac)
        sc.tl.pca(atac)

# --------------------------------------------------------------------------------------------------------------------
#                                 DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION
# --------------------------------------------------------------------------------------------------------------------

        print("\n===== NEAREST NEIGHBOR GRAPH CONSTRUCTION =====")
        print("\nConstruction of the nearest neighbor graph")
        sc.pp.neighbors(atac, n_neighbors=10, n_pcs=30)

        print("\n===== DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION=====")
        print("\nPerforming dimensionality reduction by running uniform manifold approximation and projection (UMAP)")
        sc.tl.umap(atac, spread=1.5, min_dist=.5, random_state=20)


# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

        # Visualize UMAP plot
        print("\nVisualized UMAP plot")
        sc.pl.umap(atac, color ='sample',legend_loc='on data',show=False)
        plt.savefig(os.path.join(args.results,'umap_plot_ATAC.png'))
        plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE ATAC DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
        print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
        mdata.mod['atac'] = atac
        mdata.update()
    else:
        print("ATAC modality does not exist in mdata.mod.")
    '''
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
