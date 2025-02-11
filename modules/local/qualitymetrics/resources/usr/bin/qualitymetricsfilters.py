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
import scanpy as sc                 # single-cell data processing
import matplotlib.pyplot as plt     # library for visualization
import seaborn as sns               # library for statistical data visualization
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
    This function calculates common quality control (QC) metrics for each sample and modality
    ,inspects QC plots for each computed QC metrics in each sample and filter cells and genes and ab based on QC metrics.
    """
# --------------------------------------------------------------------------------------------------------------------
#                                          LIBRARY CONFIG
# --------------------------------------------------------------------------------------------------------------------

    sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

#Define command line arguments with argparse
    parser = argparse.ArgumentParser(prog='QC_filter', usage='%(prog)s [options]', description = "QC metrics and filtering",
                                    epilog = "This function calculates common quality control (QC) metrics for each sample, inspects QC plots for each sample and filters cells based on QC plots.",
                                    )
    parser.add_argument('-ad','--input-h5mu-combined',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5mu format (including file names)")
    parser.add_argument('-d','--input-csv-doublets',metavar= 'CSV_DOUBLETS_TABLE', type=pathlib.Path, dest='input_csv_table',
                        required=True, help="paths of existing count matrix files in h5ad format (including file names)")
    parser.add_argument('-f', '--filter',dest='mt_threshold',type=float,default=15,help="parameters used to filter cells based on mithocondrial gene content")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.filtered.h5ad",
                        help="path and name of the output h5ad file")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_file = args.input_h5mu_files
    input_csv_table = args.input_csv_table
    output =args.out
    mt_threshold = args.mt_threshold

    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print(f"-File {input_h5mu_file}")
# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the MTX combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print(f"\nProcessing count matrix in folder {input_h5mu_file} ... ", end ='')

    mdata= md.read_h5ad(input_h5mu_file)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes")


# --------------------------------------------------------------------------------------------------------------------
#                                 READ DOUBLETS TABLE
# --------------------------------------------------------------------------------------------------------------------
    #print("\n===== READING DOUBLETS TABLE =====")
    input_csv_table_file=pd.read_csv(input_csv_table,index_col=0)

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
  
    gex = mdata.mod['gex']
    print(gex)
    print(mdata)

# --------------------------------------------------------------------------------------------------------------------
#                                 FILTER DOUBLETS
# --------------------------------------------------------------------------------------------------------------------
    gex.obs["doublets"] = input_csv_table_file['scDblFinder.class']
    gex.obs["doublet_score"] = input_csv_table_file['scDblFinder.score']
    #print("\n===== FILTER OUT DOUBLETS =====")
    #Filter based on doublet score
    cell_doublets =gex[gex.obs["doublets"] == 'doublet'].shape[0]
    print('''filter out {cell_doublets} cells which are doublets''')
    gex = gex[gex.obs.doublets != 'doublet', :]
    del gex.obs["doublets"]
    del gex.obs["doublet_score"]

# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    # Compute the fraction of mitochondrial, ribosomal and hemoglobin genes

    print("\n===== COMPUTE QUALITY METRICS {} =====")
    print(f"\nCompute fraction of mitochondrial, ribosomal and hemoglobin genes for {input_h5ad_file}")

    gex.var["mt"] = gex.var_names.str.startswith("MT-") | gex.var_names.str.startswith("mt-")
    gex.var["ribo"] = gex.var_names.str.startswith(("RPS", "RPL")) | gex.var_names.str.startswith(("rps", "rpl"))
    gex.var["hb"] = gex.var_names.str.contains(("^HB[^(P)]")) | gex.var_names.str.contains(("^hb[^(p)]"))

    sc.pp.calculate_qc_metrics(gex, qc_vars=["mt", "ribo", "hb"],percent_top=None, log1p=False, inplace=True)
# --------------------------------------------------------------------------------------------------------------------
#                           EVALUATE PERCENTILE
# --------------------------------------------------------------------------------------------------------------------
    percentiles = {
            'n_genes_by_counts': {
                5: round(np.percentile(gex.obs['n_genes_by_counts'], 5)),
                95: round(np.percentile(gex.obs['n_genes_by_counts'], 95))
            },
            'total_counts': {
                5: round(np.percentile(gex.obs['total_counts'], 5)),
                95: round(np.percentile(gex.obs['total_counts'], 95))
            }
        }
# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------
# Visualize quality metrics: highest expressed genes, number of genes expressed, total counts per cell and fraction of mitochondrial, ribosomal and hemoglobin genes

    fig, ax = plt.subplots(figsize=(40,10))

    print("\nVisualized the number of cells for each sample before filtering")
    sns.histplot(gex.obs, x="sample", stat="count", ax=ax)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90.)
    plt.savefig(os.path.join(args.results,'Cells_before_filtering.png'))
    plt.close()
    for sample in gex.obs['sample'].unique():
        print(f"\nVisualize density plot showing number of genes expressed, total counts per cell in {sample}")
        ax1 = plt.subplot(1, 2, 1)
        sns.histplot(gex[gex.obs['sample']== sample].obs['total_counts'], stat="count", bins=500, color='chocolate', kde=True, ax=ax1)
        plt.axvline(percentiles['total_counts'][5], color='blue', linestyle='--')
        plt.axvline(percentiles['total_counts'][95], color='blue', linestyle='--')
        ax1.set_xlim([0., 60000.])
        ax2 = plt.subplot(1, 2, 2)
        sns.histplot(gex[gex.obs['sample']== sample].obs['n_genes_by_counts'], stat="count", bins=100, color='orange', kde=True, ax=ax2)
        plt.axvline(percentiles['n_genes_by_counts'][5], color='blue', linestyle='--')
        plt.axvline(percentiles['n_genes_by_counts'][95], color='blue', linestyle='--')
        ax2.set_xlim([0., 10000.])

        plt.tight_layout()
        plt.savefig(os.path.join(args.results, f'QC_Density_{sample}.png'))
        plt.close()
        print(f"\nVisualize density plot showing fraction of mitochondrial and ribosomal genes in {sample}")
        ax1 = plt.subplot(1, 2, 1)
        sns.histplot(gex[gex.obs['sample']== sample].obs['pct_counts_mt'], stat="count", bins=100, kde=True, color='limegreen', ax=ax1)
        plt.axvline(mt_threshold, 0, 1, c='red', linestyle='--')
        ax1.set_xlim([0., 25.])
        ax2 = plt.subplot(1, 2, 2)
        sns.histplot(gex[gex.obs['sample']== sample].obs['pct_counts_ribo'], stat="count", bins=100, kde=True, color='deepskyblue', ax=ax2)
        ax2.set_xlim([0., 60.])
        plt.savefig(os.path.join(args.results, f'QC_Density_MT-Ribo_{sample}.png'))
        plt.close()
# --------------------------------------------------------------------------------------------------------------------
#                           APPLY QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    # Filter cells of low quality

    print("\n===== FILTER CELLS BASED ON QUALITY METRICS =====")
    print('Filter low quality cells on the basis of number of counts per barcode (count depth),number of genes per barcode of mitochondrial, and fraction of counts from mitochondrial genes per barcode')

    #Filter based on MIN_COUNT
    mu.pp.filter_obs(gex, 'total_counts',lambda x: x >= percentiles['total_counts'][5])
    
    #Filter based on MIN_GENES
    mu.pp.filter_obs(gex, 'n_genes_by_counts',lambda x: x >= percentiles['n_genes_by_counts'][5])

    print("Count matrix for combined samples has {} cells and {} genes after filtering".format(gex.shape[0],gex.shape[1]))
    #Filter based on MT_PERCENTAGE
    cell_number =gex[gex.obs.pct_counts_mt >= MT_PERCENTAGE].shape[0]
    print(cell_number)
    print('filter out {} cells for which the expression of mithocondrial genes is more than {}%'.format(cell_number,MT_PERCENTAGE))
    mu.pp.filter_obs(gex,'pct_counts_mt', lambda x: x < MT_PERCENTAGE)
    

    #Filter based on number of cells
    MIN_CELLS = round(1/100 * gex.shape[0])
    mu.pp.filter_var(gex, 'n_cells_by_counts', lambda x: x >= MIN_CELLS)


    print("Count matrix has {} cells and {} genes".format(gex.shape[0],gex.shape[1]))

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(20,10))
    print("\nVisualized the number of cells after filtering for each sample")
    sns.histplot(gex.obs, x="sample", stat="count", ax=ax)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90.)
    plt.savefig(os.path.join(args.results,'Cells_after_filtering.png'))
    plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
    mdata.mod['gex'] = gex
    mdata.update()
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    
    print("\n===== SAVING OUTPUT FILE =====")
    print("Saving h5ad data to file {}".format(output))
    mdata.write(output)
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
