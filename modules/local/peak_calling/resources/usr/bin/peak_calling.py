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
    This function compute the cell clutering to identify populations with similar chromatin accessibility profiles and peak calling .

    """

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

#Define command line arguments with argparse
    parser = argparse.ArgumentParser(prog='Cell clustering ', usage='%(prog)s [options]', description = "Cell clustering and peak calling for ATAC data",
                                    epilog = "This function compute the cell clutering to identify populations with similar chromatin accessibility profiles and based on clustering the peak calling.",
                                    )
    parser.add_argument('-ad','--input-h5ad-combined',metavar= 'H5AD_INPUT_FILES', type=pathlib.Path, dest='input_h5ad_files',
                        required=True, help="paths of existing matrix files in h5ad format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.peaks_atac.h5ad",
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
    print(f"MuData matrix for combined samples has {adata_atac.shape[0]} cells and {adata_atac.shape[1]} fragments")


# --------------------------------------------------------------------------------------------------------------------
#                           DIMENSIONALITY REDUCTION
# --------------------------------------------------------------------------------------------------------------------

    # Compute the spectral embedding for dimensionality reduction
    print("\n===== COMPUTE KNN GRAPH CONSTRUCTION =====")
    snap.pp.knn(adata_atac,n_neighbors=15,use_rep='X_spectral',random_state=0,inplace=True)

    # CLuster the cells using Leiden clustering
    print("\n===== CLUSTERING CELLS =====")
    print("Clustering cells using Leiden clustering ... ", end='')
    snap.tl.leiden(adata_atac, resolution=1.0, random_state=0,inplace=True)
    print("Done!")

    # Perform peak calling with MACS3
    print("\n===== COMPUTE PEAK CALLING =====")
    print("Computing peak calling ... ", end='')
    snap.tl.macs3(adata_atac,groupby='leiden',call_broad_peaks=False,inplace=True)
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           MERGE PEAKS FROM DIFFERENT GROUPS
# --------------------------------------------------------------------------------------------------------------------
    # Merge peaks from different groups
    print("\n===== MERGE PEAKS FROM DIFFERENT GROUPS =====")
    print("Merging peaks from different groups ... ", end='')
    snap.tl.merge_peaks(adata_atac.uns['macs3'], snap.genome.hg38,half_width=250)
    print("Done!")


# --------------------------------------------------------------------------------------------------------------------
#                           CREATE PEAKS COUNT MATRIX
# --------------------------------------------------------------------------------------------------------------------
 
    # Create peaks count matrix
    print("\n===== CREATE PEAKS COUNT MATRIX =====")
    print("Creating peaks count matrix ... ", end='')
    snap.pp.make_peak_matrix(adata_atac, use_rep="macs3", counting_strategy='paired-insertion', inplace=True)
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
