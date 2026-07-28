#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import pathlib                      # library for handle filesystem paths
import glob
import scanpy as sc                 # single-cell data processing
import scirpy as ir                 # single-cell AIRR-data
import anndata as ad                # store annotated matrix as anndata object
import pandas as pd
import os

warnings.filterwarnings("ignore")

# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function concatenates csv files from vdj modality.
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

    parser = argparse.ArgumentParser(prog='Concatenate_vdj', usage='%(prog)s [options]',description = "VDJ data concatenation",
                        epilog = "This function concatenated vdj filtered contig annotation files into a single csv files.")
    parser.add_argument('-ai', '--input-vdj-dir', metavar='VDJ_INPUT_FILES',nargs='+',type=pathlib.Path, dest='input_vdj_files',
                        help="paths of existing directory containing vdj matrix files in csv format (including file names)")
    parser.add_argument('-id', '--input-run-id', metavar='INPUT_RUN_ID', nargs='+', dest='input_run_id',
                        help="names of the run-id corresponding to the input adata")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="combined.vdj.h5ad",
                        help="name of the h5ad object containing the concatenated vdj table")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== VDJ FILES =====")
    input_vdj_files = args.input_vdj_files
    input_run_id = args.input_run_id
    output = args.out

    if len(input_run_id) != len(input_vdj_files):
        raise ValueError("The number of run IDs must match the number of VDJ input files.")

    # print info on the available matrices
    print("Reading vdj matrix from the following files:")
    for run, mtx in zip(input_run_id, input_vdj_files):
        print(f"Run: {run:15s} - File: {mtx}")


# --------------------------------------------------------------------------------------------------------------------
#                                 READ VDJ FILES
# --------------------------------------------------------------------------------------------------------------------

    adata_vdj_list = []

    for run, vdj_path in zip(input_run_id, input_vdj_files):
        vdj = pathlib.Path(vdj_path)
        print("\n===== READING CONTIGUE ANNOTATION MATRIX =====")
        print("\nProcessing filtered contigue table in file ... ", end ='')
        adata_vdj = ir.io.read_10x_vdj(str(vdj), filtered=False)
        n_cells = adata_vdj.n_obs
        print(f"Done! Number of cells for sample {run}: {n_cells}")
        adata_vdj.obs['sample'] = run
        adata_vdj_list.append(adata_vdj)

    if not adata_vdj_list:
        print("No valid input file provided. Skipping reading of the vdj annotation.")

# --------------------------------------------------------------------------------------------------------------------
#                           VDJ TABLE CONCATENATION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== CONCATENATING VDJ TABLES =====")

    if len(adata_vdj_list) == 0:
        print("No valid files were found. Nothing to save.")
    else:
        if len(adata_vdj_list) == 1:
            adata_vdj_concatenated = adata_vdj_list[0]
            print("Only one non-empty file found. Saving the file as is without concatenation.")
        else:
            adata_vdj_concatenated = ad.concat(
                adata_vdj_list,
                join="outer",
                merge="same",
                index_unique="_"
            )

        # Ensure sample labels are correctly populated after concatenation.
        if "sample" not in adata_vdj_concatenated.obs.columns:
            sample_labels = []
            for adata_vdj in adata_vdj_list:
                run_label = adata_vdj.obs["sample"].iloc[0] if "sample" in adata_vdj.obs.columns else None
                sample_labels.extend([run_label] * adata_vdj.n_obs)
            adata_vdj_concatenated.obs["sample"] = sample_labels

        print(f"Concatenated vdj table for {len(adata_vdj_list)} files has {adata_vdj_concatenated.shape[0]} cells")
        print("Done!")
        print("Number of cells per sample after concatenation:")
        counts_per_sample = adata_vdj_concatenated.obs['sample'].value_counts()
        print(counts_per_sample)


        print("\n===== SAVING OUTPUT FILE =====")

        print(f"Saving vdj table data in {output}")
        adata_vdj_concatenated.write(output)
        print("Done!")

#####################################################################################################


if __name__ == '__main__':
    main()
