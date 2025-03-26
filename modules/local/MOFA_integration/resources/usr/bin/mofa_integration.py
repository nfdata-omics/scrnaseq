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
    output = args.out


    # print info on the available matrix
    print("Reading combined count matrix from the following file:")
    print(f"-File {str(input_h5mu_files)}:")

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
#                                 MOFA CALCULATION ON SELECTED MODALITIES
# --------------------------------------------------------------------------------------------------------------------
    gex= mdata.mod['gex']
    pro= mdata.mod['pro']
    modalities = {}
    modalities["gex"] = gex
    modalities["pro"] = pro
    mdata_subset = MuData(modalities)

    print("\n===== MOFA CALCULATION =====")
    # MOFA calculation
    #sc.pp.neighbors(mdata_subset['gex'])
    #sc.pp.neighbors(mdata_subset['pro'])
    print("\nCalculating MOFA ... ")
    mu.tl.mofa(mdata_subset,use_obs='union')
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
    mu.pl.embedding(mdata_subset, basis="umap_mofa",color ='sample',legend_loc="on data", show=False)
    plt.savefig(os.path.join(args.results,'umap_plot_mofa.png'))
    plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
    mdata.mod['gex'] = mdata_subset.mod['gex']
    mdata.mod['pro'] = mdata_subset.mod['pro']
    mdata.update()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    print(f"Saving h5mu data to file {output}")
    mdata_subset.write(output)
    print("Done!")


#####################################################################################################


if __name__ == '__main__':
    main()
    