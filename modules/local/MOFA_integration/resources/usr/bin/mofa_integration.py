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
import scanpy as sc                 # single-cell data processing
import muon as mu
from mudata import MuData
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
    This function computed the cell neigbourhood graph for all modalities using a weighted nearest neighbours (WNN) method.
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

    parser = argparse.ArgumentParser(prog='LogNorm', usage='%(prog)s [options]',description = "Multi-Omics Factor Analysis (MOFA) integration of count matrix for each modality",
                        epilog = "This function learns an interpretable latent space jointly on multiple modalities")
    parser.add_argument('-ad','--input-h5mu-file',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing matrix files in h5mu format (including file names)")
    parser.add_argument('-at', '--input-atac-file', metavar='ATAC_INPUT_FILES',type=pathlib.Path, dest='input_atac_files',
                        default=pathlib.Path(''),help="paths of existing atac matrix files in h5ad format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.mofa.h5mu",
                        help="path and name of the output h5mu file after filtering")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()
    parser.print_help()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_files = args.input_h5mu_files
    input_atac_file = args.input_atac_files
    output = args.out


    # print info on the available matrix
    print("Reading combined count matrix from the following file:")
    print(f"-File {str(input_h5mu_files)}:")

    print("Reading combined atac count matrix from the following file:")
    print(f"-File {input_atac_file}")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5MU FILES
# --------------------------------------------------------------------------------------------------------------------


     # Read folders with the combined count matrice and store datasets in a dictionary

    print("\n===== READING COMBINED H5MU MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing filtered count matrix in folder ... ", end ='')
    mdata= md.read(input_h5mu_files)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")


# --------------------------------------------------------------------------------------------------------------------
#                                 READ ATAC FILES
# --------------------------------------------------------------------------------------------------------------------
    if input_atac_file and input_atac_file != pathlib.Path(''):
        # Read folders with the ATAC combined count matrice and store datasets in a dictionary
        print("\n===== READING COMBINED ATAC MATRIX =====")
        # read the gex count matrix for the combined samples and print some initial info
        print("\nProcessing count matrix in folder ... ", end ='')
        adata_atac= sc.read_h5ad(input_atac_file)
        print("Done!")
        print(f"Atac count matrix for combined samples has {adata_atac.shape[0]} cells and {adata_atac.shape[1]} peaks")

# --------------------------------------------------------------------------------------------------------------------
#                                 ADD ATAC MODALITY TO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------

        mdata.mod['atac'] = adata_atac
        print("ATAC modality added to MuData")

# --------------------------------------------------------------------------------------------------------------------
#                                 MOFA CALCULATION ON SELECTED MODALITIES
# --------------------------------------------------------------------------------------------------------------------
    modalities = {}
    modalities["gex"] = mdata.mod['gex']
    if 'pro' in mdata.mod:
        modalities["pro"] = mdata.mod['pro']
    if 'atac' in mdata.mod:
        modalities["atac"] = mdata.mod['atac']
    mdata_subset = MuData(modalities)


    for col in mdata_subset.var.select_dtypes(include="boolean"):
        mdata_subset.var[col] = mdata_subset.var[col].fillna(False).astype(bool)

    for col in mdata_subset.obs.select_dtypes(include="boolean"):
        mdata_subset.obs[col] = mdata_subset.obs[col].fillna(False).astype(bool)




    print("\n===== MOFA CALCULATION =====")
    # MOFA calculation
    print("\nCalculating MOFA ... ")
    mu.tl.mofa(mdata_subset,use_obs='union',seed=42,use_var='highly_variable')
    sc.pp.neighbors(mdata_subset, use_rep="X_mofa")
    print("Done!")

    print("\n===== DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION=====")
    print("\nPerforming dimensionality reduction by running uniform manifold approximation and projection (UMAP)")
    sc.tl.umap(mdata_subset)
    mdata_subset.obsm["X_umap_mofa"] = mdata_subset.obsm["X_umap"].copy()
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(20,10))
    # Visualize UMAP plot
    print("\nVisualized UMAP plot")
    mu.pl.embedding(mdata_subset, basis="umap_mofa",color = 'sample',legend_loc="on data", show=False)
    plt.savefig(os.path.join(args.results,'umap_mofa.png'))
    plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
    mdata.mod['gex'] = mdata_subset.mod['gex']
    if 'pro' in mdata_subset.mod:
        mdata.mod['pro'] = mdata_subset.mod['pro']
    if 'atac' in mdata_subset.mod:
        mdata.mod['atac'] = mdata_subset.mod['atac']
    mdata.update()

    cols_to_remove = ["mt", "ribo", "hb", "gene_pass", "highly_variable"]
    mdata.var = mdata.var.drop(columns=[c for c in cols_to_remove if c in mdata.var.columns])

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    print(f"Saving h5mu data to file {output}")
    mdata.write(output)
    print("Done!")


#####################################################################################################


if __name__ == '__main__':
    main()
