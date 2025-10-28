#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import warnings
import argparse                     # command line arguments parser
import pathlib                      # library for handle filesystem paths
import os
import scanpy as sc                 # single-cell data processing
import mudata as md
import muon as mu
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function lognormalized the matrices for each modality.
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

    parser = argparse.ArgumentParser(prog='LogNorm', usage='%(prog)s [options]',description = "Normalization and logaritmic trasformation of count matrix for each modality",
                        epilog = "This function normalize and logarithmize the data for each modality")
    parser.add_argument('-ad','--input-h5mu-file',metavar= 'H5MU_INPUT_FILES',type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing matrix files in h5mu format (including file names)")
    parser.add_argument('-rw','--input-h5ad-raw-file',metavar= 'H5AD_RAW_INPUT_FILES', type=pathlib.Path, dest='input_h5ad_raw_files',
                        default=pathlib.Path(""),help="paths of existing raw matrix files in h5ad format (including file names)")
    parser.add_argument('-cycle', '--input_cellcycle_file',metavar='CELLCYCLE_FILE', type=pathlib.Path, dest='cellcycle_file',
                    help="Path of an csv input file (.csv)")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.norm.h5mu",
                        help="path and name of the output h5mu file after filtering")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()
    parser.print_help()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_files = args.input_h5mu_files
    input_h5ad_files = args.input_h5ad_raw_files
    output = args.out
    input_cellcycle = args.cellcycle_file


    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print(f"-File {str(input_h5mu_files)}:")
    if input_h5ad_files and input_h5ad_files != pathlib.Path(""):
        print("Reading raw count matrix from the following file:")
        print(f"-File {str(input_h5ad_files)}:")
    else:
        print("No raw H5AD file provided, skipping raw/CITE processing.")


# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5MU AND H5AD FILES
# --------------------------------------------------------------------------------------------------------------------


    # Read folders with the combined count matrice and store datasets in a dictionary

    print("\n===== READING COMBINED H5MU MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing filtered count matrix in folder ... ", end ='')

    mdata= md.read(input_h5mu_files)

    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD RAW FILES
# --------------------------------------------------------------------------------------------------------------------
    # Read folders with the combined raw count matrice and store datasets in a dictionary

    print("\n===== READING COMBINED H5AD MATRIX =====")
    # read the raw count matrix for the combined samples and print some initial info
    print("\nProcessing raw count matrix in folder ... ", end ='')

    if input_h5ad_files and input_h5ad_files != pathlib.Path(""):
        adata_raw = sc.read_h5ad(input_h5ad_files)
        if 'feature_types' in adata_raw.var:
            pro_raw = adata_raw[:, adata_raw.var["feature_types"] == "Antibody Capture"].copy()
            print("Done!")
            print(f"Raw count matrix for combined samples has {pro_raw.shape[0]} cells and {pro_raw.shape[1]} genes/ab")
        else:
            print("No Antibody Capture features found in raw H5AD.")
            pro_raw = None
    else:
        adata_raw = None
        pro_raw = None
        print("No raw H5AD file provided.")

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== GEX MODALITY DATA =====")
    gex = mdata.mod['gex']

# --------------------------------------------------------------------------------------------------------------------
#                                 NORMALIZATION
# --------------------------------------------------------------------------------------------------------------------
    # Saving count data before normalization
    print("Saving count data before normalization in slot Count.")
    gex.layers["count"] = gex.X.copy()

    print("\n===== NORMALIZATION =====")
    # Normalizing to median total counts
    print("\nNormalize to median total counts ... ")
    sc.pp.normalize_total(gex)

    print("Done!")

    print("\n===== LOGARITMIC TRASFORMATION =====")
    print("\nLogarithmize the data ... ", end ='')
    # Logarithmize the data
    sc.pp.log1p(gex)

    print("Done!")
    gex.layers["normalized_gex"] = gex.X.copy()

# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE CELL CYCLE SCORE
# --------------------------------------------------------------------------------------------------------------------
    print(f"Reading cell cycle gene list from: {input_cellcycle}")

    with open(input_cellcycle) as f:
        cell_cycle_genes = [x.strip() for x in f]

    original_var_names = gex.var_names.copy()
    gex.var['gene_name'] = gex.var['gene_name'].astype(str)
    gex.var_names = gex.var['gene_name']
    gex.var_names_make_unique()
    # Define genes associated to the S phase and genes associated to the G2M phase
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in gex.var_names]
    print("\n===== COMPUTING CELL CYCLE SCORE =====")
    print("Computing cell cycle score ... ", end='')
    sc.tl.score_genes_cell_cycle(gex, s_genes=s_genes, g2m_genes=g2m_genes)
    print("Done!")

    #Compute pca on subset of genes
    print("\nComputing PCA on cell cycle genes ... ", end='')
    gex_cc_genes = gex[:, cell_cycle_genes]
    sc.tl.pca(gex_cc_genes)
    sc.pl.pca_scatter(gex_cc_genes, color='phase',show=False)
    plt.savefig(os.path.join(args.results,'pca_cellcycle_GEX_phase.pdf'), bbox_inches='tight', dpi=300)
    print("Done!")
    plt.close()

    # Visualize PCA plot
    print("\nVisualized PCA plot")
    plt.figure(figsize=(14, 13))
    sc.pl.pca(gex_cc_genes, color='sample', show=False)
    plt.savefig(os.path.join(args.results,'pca_cellcycle_GEX_sample.pdf'), bbox_inches='tight', dpi=300)
    plt.close()

    # Add cell cycle score to the gex object
    gex.obs['S_score'] = gex.obs['S_score'].astype(float)
    gex.obs['G2M_score'] = gex.obs['G2M_score'].astype(float)

    # Save gex object into mdata
    gex.var_names = original_var_names


    mdata.mod['gex'] = gex
    mdata.update()

# --------------------------------------------------------------------------------------------------------------------
#                                 CITE MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    if 'pro' in mdata.mod:
        print("\n===== CITE MODALITY DATA =====")
        pro = mdata.mod['pro']

# --------------------------------------------------------------------------------------------------------------------
#                                 NORMALIZATION
# --------------------------------------------------------------------------------------------------------------------
        # Saving count data before normalization
        print("Saving count data before normalization in slot Count.")
        pro.layers["count"] = pro.X.copy()
        print("\n===== NORMALIZATION =====")
        # Denoising and normalizing protein expression with DSB (Denoised and Scaled by Background)
        print("\nDenoising and normalize with Denoised and Scaled by Background method... ")
        mu.prot.pp.dsb(pro, pro_raw)

        pro.layers["normalized_pro"] = pro.X.copy()
        mdata.mod['pro'] = pro
        mdata.update()

        print("Done!")
    else:
        print("CITE modality does not exist in mdata.mod.")


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
