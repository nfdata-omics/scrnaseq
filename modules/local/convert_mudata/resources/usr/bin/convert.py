#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import pathlib                      # library for handle filesystem paths
import scanpy as sc                 # single-cell data processing
import pandas as pd
from mudata import MuData


warnings.filterwarnings("ignore")

# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function creates a MuData object.
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

    parser = argparse.ArgumentParser(prog='Create MuData object', usage='%(prog)s [options]',description = "MuData object convertion",
                        epilog = "This function creates a MuData object for storing GEX,VDJ and CITE-seq data.")
    parser.add_argument('-ad','--input-file',metavar= 'INPUT_FILES', type=pathlib.Path, dest='input_files',
                        help="paths of existing count matrix files in h5ad format (including file names)")
    parser.add_argument('-ai', '--input-vdj-file', metavar='VDJ_INPUT_FILES',type=pathlib.Path, dest='input_vdj_files',
                        default=pathlib.Path(''),help="paths of existing vdj matrix files in h5ad format (including file names)")
    parser.add_argument('-csv','--input-csv-file',metavar= 'CSV_INPUT_FILES', type=pathlib.Path,  dest='input_csv_files',
                        help="paths of existing metadata table in csv format")
    parser.add_argument('-meta','--metadata-file',metavar= 'METADATA_INPUT_FILES', type=pathlib.Path,  dest='input_metadata_files',
                        help="paths of existing metadata table in csv format")
    parser.add_argument('-o', '--out', metavar='MUDATA_OUTPUT_FILE', type=pathlib.Path, default="matrix.mudata.h5mu",
                        help="name of the muData object")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT GEX and VDJ FILES =====")
    input_file = args.input_files
    input_vdj_file = args.input_vdj_files
    input_csv_file = args.input_csv_files
    input_metadata_file = args.input_metadata_files
    output = args.out

    # print info on the available matrices
    print("Reading combined gex count matrix from the following file:")
    print(f"-File {input_file}")

    print("Reading filtered annotation table from the following file:")
    print(f"-File {input_vdj_file}")

    print("Reading metadata table from the following file:")
    print(f"-File {input_csv_file}")

    print("Reading sample metadata information from the following file:")
    print(f"-File {input_metadata_file}")



# --------------------------------------------------------------------------------------------------------------------
#                                 READ GEX AND AB FILES
# --------------------------------------------------------------------------------------------------------------------
    if input_file:
        # Read folders with the MTX combined count matrice and store datasets in a dictionary
        print("\n===== READING COMBINED MATRIX =====")
        # read the gex count matrix for the combined samples and print some initial info
        print("\nProcessing count matrix in folder ... ", end ='')
        adata= sc.read_h5ad(input_file)
        print("Done!")
        print(f"Gex count matrix for combined samples has {adata.shape[0]} cells and {adata.shape[1]} genes")
    else:
        print("No valid input file provided. Skipping reading of the count matrix.")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ VDJ FILES
# --------------------------------------------------------------------------------------------------------------------
    if input_vdj_file and input_vdj_file != pathlib.Path(''):
        # Read folders with the filtered contigue annotation and store datasets in a dictionary
        print("\n===== READING CONTIGUE ANNOTATION MATRIX =====")
        # read the filtered contigue annotation file for the combined samples and print some initial info
        print("\nProcessing filtered contigue table in folder ... ", end ='')
        adata_vdj= sc.read_h5ad(input_vdj_file)
        print("Done!")
    else:
        print("No valid input file provided. Skipping reading of the vdj annotation.")



    adata_vdj = None


# --------------------------------------------------------------------------------------------------------------------
#                                 READ CSV FILES
# --------------------------------------------------------------------------------------------------------------------
    metadata_df = None
    if input_csv_file and input_csv_file.exists():
        print("\n===== READING METADATA CSV =====")
        metadata_df = pd.read_csv(input_csv_file,sep='\t',header=0)
        metadata_df['pool_barcode'] = metadata_df['Barcode'].astype(str) + '_' + metadata_df['pool'].astype(str) + '_cellbender_filter'
        print(metadata_df)
        print(f"Metadata table has {metadata_df.shape[0]} rows and {metadata_df.shape[1]} columns")
    else:
        print("No valid CSV file provided. Skipping reading of the metadata table.")

    meta_df = None
    if input_metadata_file and input_metadata_file.exists() and input_metadata_file.stat().st_size > 0:
        print("\n===== READING SAMPLE METADATA CSV =====")
        meta_df = pd.read_csv(input_metadata_file, sep=',', header=0)
        print(meta_df)
        print(f"Sample metadata table has {meta_df.shape[0]} rows and {meta_df.shape[1]} columns")
    else:
        print("No valid metadata CSV file provided. Skipping reading of the sample metadata table.")

# --------------------------------------------------------------------------------------------------------------------
#                                 ADDED METADATA TO OBS
# --------------------------------------------------------------------------------------------------------------------

    if metadata_df is not None:
        print("\n===== ADDING METADATA TO OBS =====")

        if 'Barcode' in metadata_df.columns:
            metadata_df = metadata_df.set_index('pool_barcode')
            obs_names_index = pd.Index(adata.obs_names)
            intersect_barcodes = obs_names_index.intersection(metadata_df.index)
            metadata_to_add = metadata_df.loc[intersect_barcodes]
            adata.obs = adata.obs.join(metadata_to_add, how='left')
            adata.obs['Souporcell_Cluster'] = adata.obs['Souporcell_Cluster'].astype(str)
            print(f"Metadata joined to MuData obs for {len(intersect_barcodes)} barcodes.")
        else:
            print("No 'Barcode' column found in metadata CSV; skipping join.")

    if meta_df is not None:
        print("\n===== ADDING SAMPLE METADATA TO OBS =====")

        if 'sample' in meta_df.columns:
            adata.obs['sample'] = adata.obs['sample'].astype(str).str.replace('_filtered', '', regex=False).str.replace('_parse', '', regex=False).str.strip()
            meta_df['sample'] = meta_df['sample'].astype(str)
            adata.obs = adata.obs.join(meta_df.set_index('sample'), on='sample', how='left').astype(str)
            print("Sample metadata joined to MuData obs.")
        else:
            print("No 'sample' column found in metadata CSV; skipping join.")
# --------------------------------------------------------------------------------------------------------------------
#                                 CREATE MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    #Creates dictionary to store all modalities
    modalities = {}
    try:
        # Add 'gex' modality if defined
        if adata[:, adata.var["feature_types"] == "Gene Expression"].shape[1] > 0:
            modalities["gex"] = adata[:, adata.var["feature_types"] == "Gene Expression"]
        # Add 'pro' modality if defined
        if adata[:, adata.var["feature_types"] == "Antibody Capture"].shape[1] > 0:
            modalities["pro"] = adata[:, adata.var["feature_types"] == "Antibody Capture"]
    except NameError:
        pass

    try:
        # Add 'airr' modality if defined
        if adata_vdj is not None:
            modalities["airr"] = adata_vdj
    except NameError:
        pass

    # Creates MuData object
    mdata = MuData(modalities)
    if 'airr' in mdata.mod:
        mdata.obs['airr:sample'] = mdata.obs['airr:sample'].astype(str)
    mdata.update()


# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== SAVING OUTPUT FILE =====")

    print(f"Saving MuData object to file {output}")
    mdata.write(output)
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
