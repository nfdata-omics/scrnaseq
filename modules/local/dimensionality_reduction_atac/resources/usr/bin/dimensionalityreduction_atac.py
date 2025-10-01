#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import argparse                     # command line arguments parser
import warnings
import os                           # filesystem utilities
import pathlib                      # library for handle filesystem paths
import numpy as np
import pandas as pd                 # library for data analysis and manipulation                
import snapatac2 as snap
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
    This function compute the dimensionality reduction for ATAC data using TF-IDF normalization and PCA
    """

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

#Define command line arguments with argparse
    parser = argparse.ArgumentParser(prog='DimensionalityRed', usage='%(prog)s [options]', description = "Dimensionality reduction for ATAC data",
                                    epilog = "This function compute the dimensionality reduction using TF-IDF normalization.",
                                    )
    parser.add_argument('-ad','--input-h5ad-combined',metavar= 'H5AD_INPUT_FILES', type=pathlib.Path, dest='input_h5ad_files',
                        required=True, help="paths of existing matrix files in h5ad format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.dimred_atac.h5ad",
                        help="path and name of the output h5ad file")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5ad_file = args.input_h5ad_files
    output =args.out
    
    # print info on the available matrices
    print("Reading fragment file from the following file:")
    print(f"-File {input_h5ad_file}")

    
# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print(f"\nProcessing AnnData object in folder {input_h5ad_file} ... ", end ='')

    adata_atac= snap.read(input_h5ad_file,backed=None)
    print(adata_atac.obs)
    print(adata_atac.var)
    print("Done!")
    print(f"MuData matrix for combined samples has {adata_atac.shape[0]} cells and {adata_atac.shape[1]} fragments")


# --------------------------------------------------------------------------------------------------------------------
#                           DIMENSIONALITY REDUCTION
# --------------------------------------------------------------------------------------------------------------------

    # Compute the spectral embedding for dimensionality reduction
    print("\n===== COMPUTE SPECTRAL EMBEDDING =====")
    snap.pp.add_tile_matrix(adata_atac,bin_size=500,exclude_chroms=["chrM"],min_frag_size=50,max_frag_size=1000,counting_strategy='paired-insertion',inplace=True)


    # Normalize the data using TF-IDF normalization
    print("\nNormalizing data using TF-IDF normalization ... ", end='')
    snap.pp.select_features(adata_atac, n_features=500000)
    print(adata_atac)
    snap.tl.spectral(adata_atac,n_comps=30,features="selected",random_state=0,inplace=True)
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           BATCH CORRECTION
# --------------------------------------------------------------------------------------------------------------------

    # Perform batch correction
    # Importantly set groupby='variable to preserve' parameter if you want to preserve that difference in the batch correction
    print("\n===== PERFORM BATCH CORRECTION =====")
    print("Performing batch correction ... ", end='')
    snap.pp.mnc_correct(adata_atac, batch="sample", key_added='X_spectral')
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           CLUSTERING ON TILE MATRIX
# --------------------------------------------------------------------------------------------------------------------

    # Perform clustering
    print("\n===== PERFORM CLUSTERING =====")
    print("Performing clustering ... ", end='')
    snap.pp.knn(adata_atac)
    snap.tl.leiden(adata_atac,key_added='leiden_tile',inplace=True) 
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE AND VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------
    
    # Compute UMAP for visualization
    print("\n===== COMPUTE UMAP =====")
    print("Computing UMAP ... ", end='')
    #snap.tl.umap(adata_atac)
    print("Done!")

    # Visualize UMAP plot
    print("\nVisualized UMAP plot")
    #snap.pl.umap(adata_atac, color="sample",interactive=False)
    print("Done!")
    
    
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    print(f"Saving atac data in {output}")
    adata_atac.write(output)
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
