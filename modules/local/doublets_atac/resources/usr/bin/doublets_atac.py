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




warnings.filterwarnings("ignore")
# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function calculates common quality control (QC) metrics for each sample and ATAC modality
    ,inspects QC plots for each computed QC metrics in each sample and filter cells based on QC metrics.
    """

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

#Define command line arguments with argparse
    parser = argparse.ArgumentParser(prog='QC_filter', usage='%(prog)s [options]', description = "QC metrics and filtering",
                                    epilog = "This function calculates common quality control (QC) metrics for each sample and ATAC modality, inspects QC plots for each sample and filters cells based on QC plots.",
                                    )
    parser.add_argument('-ad','--input-h5ad-combined',metavar= 'H5AD_INPUT_FILES', type=pathlib.Path, dest='input_h5ad_files',
                        required=True, help="paths of existing matrix files in h5ad format (including file names)")
    #parser.add_argument('-id', '--input-run-id', metavar='INPUT_RUN_ID', dest='input_run_id',
    #                    help="names of the run-id corresponding to the input adata")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.doublets_atac.h5ad",
                        help="path and name of the output h5ad file")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5ad_file = args.input_h5ad_files
    #input_run_id = args.input_run_id
    output =args.out
    
    # print info on the available matrices
    print("Reading fragment file from the following file:")
    print(f"-File {input_h5ad_file}")

    #for run, fragment in zip(input_run_id, input_fragment_file):
    #print(f"Run: {run:15s} - File: {fragment}")
    
# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print(f"\nProcessing AnnData object in folder {input_h5ad_file} ... ", end ='')

    adata_atac= snap.read_dataset(input_h5ad_file)
    print("Done!")
    print(f"MuData matrix for combined samples has {adata_atac.shape[0]} cells and {adata_atac.shape[1]} fragments")

# --------------------------------------------------------------------------------------------------------------------
#                           CELL BY BIN MATRIX
# --------------------------------------------------------------------------------------------------------------------

    # Compute the cell by bin matrix for each sample
    print("\n===== COMPUTE CELL BY BIN MATRIX =====")
    snap.pp.add_tile_matrix(adata_atac,bin_size=500)
    snap.pp.select_features(adata_atac, n_features=100000)
    snap.pp.scrublet(adata_atac)
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
