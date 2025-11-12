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
    parser.add_argument('-meta', '--input-metadata',metavar='META_FILE', type=pathlib.Path, dest='meta_file',
                    help="Path of an csv input file (.csv)")
    parser.add_argument('-t', '--tile_out', metavar='H5AD_OUTPUT_TILE', type=pathlib.Path, default="matrix.tile_atac.h5ad",
                        help="path and name of the output h5ad file with the tile matrix")
    parser.add_argument('-p', '--peak_out', metavar='H5AD_OUTPUT_PEAK', type=pathlib.Path, default="matrix.peak_atac.h5ad",
                        help="path and name of the output h5ad file with the peak matrix")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5ad_file = args.input_h5ad_files
    input_meta_file = args.meta_file
    output_tile = args.tile_out
    output_peak = args.peak_out

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
    print(f"MuData matrix for combined samples has {adata_atac.shape[0]} cells and {adata_atac.shape[1]} fragments")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ META DATA
# --------------------------------------------------------------------------------------------------------------------
    
    # # Read metadata from an Excel file if provided
    # if input_meta_file:
    #     print("\n===== READING META DATA =====")
    #     print(f"Reading metadata from {input_meta_file} ... ", end='')
    #     meta_data = pd.read_csv(input_meta_file)
    #     print(meta_data)
    #     print("Done!")


    #     def transform_id(x):
    #         # remove '_filtered'
    #         x = x.replace('_filtered', '')
    #         # split su '_'
    #         parts = x.split('_')
    #         # rebuild ad  "sample:cellID"
    #         # parts[-1] is the sample, the remaining is the cellID
    #         sample = parts[-1]
    #         cell_id = "_".join(parts[:-1])
    #         return f"{sample}:{cell_id}"


    #     meta_data['cell_id'] = meta_data['Unnamed: 0'].apply(transform_id)


    #     meta_data = meta_data.set_index('cell_id')

    #     print(meta_data.head())

    #     meta_data.index = meta_data.index.astype(str)
    #     adata_atac.obs.index = adata_atac.obs.index.astype(str)

    #     # Find the common cells between adata_atac and metadata
    #     common_cells = adata_atac.obs.index.intersection(meta_data.index)
    #     print(f"\nFound {len(common_cells)} common cells between AnnData and metadata.")

    #     # Keep only the common cells in adata_atac
    #     adata_atac = adata_atac[common_cells].copy()

    #     # Filter metadata to include only those common cells
    #     meta_data_filtered = meta_data.loc[meta_data.index.isin(common_cells)]

    #     # Join metadata into adata_atac.obs
    #     adata_atac.obs = adata_atac.obs.join(meta_data_filtered, how='left')

    #     print("Done!")

    # else:
    #     print("No metadata file provided. Continuing without metadata.")


# --------------------------------------------------------------------------------------------------------------------
#                           PEAK CALLING
# --------------------------------------------------------------------------------------------------------------------

    # Perform peak calling with MACS3
    print("\n===== COMPUTE PEAK CALLING =====")
    print("Computing peak calling ... ", end='')
    # adata_atac.obs['gex:celltypist:Immune_All_High:majority_voting'] = (
    # adata_atac.obs['gex:celltypist:Immune_All_High:majority_voting']
    # .astype(str)
    # .str.replace('/', '_')
    # )
    snap.tl.macs3(adata_atac,groupby='leiden_tile',call_broad_peaks=False,inplace=True)
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           MERGE PEAKS FROM DIFFERENT GROUPS
# --------------------------------------------------------------------------------------------------------------------
    # Merge peaks from different groups
    print("\n===== MERGE PEAKS FROM DIFFERENT GROUPS =====")
    print("Merging peaks from different groups ... ", end='')
    merged_peaks = snap.tl.merge_peaks(adata_atac.uns['macs3'], snap.genome.hg38,half_width=250)
    print("Done!")

    peaks_list = merged_peaks['Peaks'].to_list()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE TILE MATRIX AS LAYER
# --------------------------------------------------------------------------------------------------------------------

    #print("\n===== SAVE TILE MATRIX AS LAYER =====")
    #print("Saving tile matrix as layer ... ", end='')
    #adata_atac.layers["tile"] = adata_atac.X.copy()
    #print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           CREATE PEAKS COUNT MATRIX
# --------------------------------------------------------------------------------------------------------------------

    # Create peaks count matrix
    print("\n===== CREATE PEAKS COUNT MATRIX =====")
    print("Creating peaks count matrix ... ", end='')
    peak_matrix = snap.pp.make_peak_matrix(adata_atac,use_rep=peaks_list,inplace=False)
    print("Done!")
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    print(f"Saving peak matrix in {output_peak}")
    peak_matrix.write(output_peak)
    print("Done!")

    print(f"Saving tile matrix in {output_tile}")
    adata_atac.write(output_tile)
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
