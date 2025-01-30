#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import argparse                     # command line arguments parser
import os                           # filesystem utilities
import scanpy as sc                 # single-cell data processing
import anndata as ad                # store annotated matrix as anndata object
import matplotlib.pyplot as plt     # library for visualization
import seaborn as sns               # library for statistical data visualization
import pandas as pd                 # library for data analysis and manipulation 
import pathlib                      # library for handle filesystem paths

# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function calculates common quality control (QC) metrics for each sample 
    ,inspects QC plots for each computed QC metrics in each sample and filter cells and genes based on QC metrics.
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

    parser = argparse.ArgumentParser(prog='QC_filter', usage='%(prog)s [options]', description = "QC metrics and filtering",
                                    epilog = "This function calculates common quality control (QC) metrics for each sample, inspects QC plots for each sample and filters cells based on QC plots.",
                                    )
    parser.add_argument('-ad','--input-h5ad-combined',metavar= 'H5AD_INPUT_FILES', type=pathlib.Path, dest='input_h5ad_files',
                        required=True, help="paths of existing count matrix files in h5ad format (including file names)")
    parser.add_argument('-d','--input-csv-doublets',metavar= 'CSV_DOUBLETS_TABLE', type=pathlib.Path, dest='input_csv_table',
                        required=True, help="paths of existing count matrix files in h5ad format (including file names)")                    
    parser.add_argument('-f', '--filter',dest='MT_PERCENTAGE',type=float,default=15,help="parameters used to filter cells based on mithocondrial gene content")
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
    input_h5ad_file = args.input_h5ad_files
    input_csv_table = args.input_csv_table
    output =args.out
    MT_PERCENTAGE = args.MT_PERCENTAGE

    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print("-File {}".format(str(input_h5ad_file)))


# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES 
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the MTX combined count matrice and store datasets in a dictionary
    
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder {} ... ".format(input_h5ad_file), end ='')

    adata= sc.read_h5ad(input_h5ad_file)
        
    print("Done!")
    print("Count matrix for combined samples has {} cells and {} genes".format(adata.shape[0],adata.shape[1]))


# --------------------------------------------------------------------------------------------------------------------
#                                 READ DOUBLETS TABLE
# --------------------------------------------------------------------------------------------------------------------
    
    #print("\n===== READING DOUBLETS TABLE =====")
    input_csv_table_file=pd.read_csv(input_csv_table,index_col=0)
    

# --------------------------------------------------------------------------------------------------------------------
#                                 FILTER DOUBLETS
# --------------------------------------------------------------------------------------------------------------------
    adata.obs["doublets"] = input_csv_table_file['scDblFinder.class']
    adata.obs["doublet_score"] = input_csv_table_file['scDblFinder.score']
    
    #print("\n===== FILTER OUT DOUBLETS =====")
    #Filter based on doublet score
    cell_doublets =adata[adata.obs["doublets"] == 'doublet'].shape[0]
    print('''filter out {} cells which are doublets'''.format(cell_doublets))
    adata = adata[adata.obs.doublets != 'doublet', :]

# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    # Compute the fraction of mitochondrial, ribosomal and hemoglobin genes 

    print("\n===== COMPUTE QUALITY METRICS {} =====")
    print("\nCompute fraction of mitochondrial, ribosomal and hemoglobin genes for {}".format(input_h5ad_file))

    adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-") 
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL")) | adata.var_names.str.startswith(("rps", "rpl"))
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]")) | adata.var_names.str.contains(("^hb[^(p)]"))  

    
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"],percent_top=None, log1p=False, inplace=True)
    
# --------------------------------------------------------------------------------------------------------------------
#                           EVALUATE PERCENTILE
# --------------------------------------------------------------------------------------------------------------------
    percentiles = {
            'n_genes_by_counts': {
                5: round(np.percentile(adata.obs['n_genes_by_counts'], 5)),
                95: round(np.percentile(adata.obs['n_genes_by_counts'], 95))
            },
            'total_counts': {
                5: round(np.percentile(adata.obs['total_counts'], 5)),
                95: round(np.percentile(adata.obs['total_counts'], 95))
            }
        }
# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    
# Visualize quality metrics: highest expressed genes, number of genes expressed, total counts per cell and fraction of mitochondrial, ribosomal and hemoglobin genes

    fig, ax = plt.subplots(figsize=(40,10))

    print("\nVisualized the number of cells for each sample before filtering")
    sns.histplot(adata.obs, x="sample", stat="count", ax=ax)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90.)
    plt.savefig(os.path.join(args.results,'Cells_before_filtering.png'))
    plt.close() 
        
    for sample in adata.obs['sample'].unique():
            
            print("\nVisualize density plot showing number of genes expressed, total counts per cell in {}".format(sample))
            ax1 = plt.subplot(1, 2, 1)
            sns.histplot(adata[adata.obs['sample']== sample].obs['total_counts'], stat="count", bins=500, color='chocolate', kde=True, ax=ax1)
            plt.axvline(percentiles['total_counts'][5], color='blue', linestyle='--')
            plt.axvline(percentiles['total_counts'][95], color='blue', linestyle='--')
            ax1.set_xlim([0., 60000.])
         
            ax2 = plt.subplot(1, 2, 2)
            sns.histplot(adata[adata.obs['sample']== sample].obs['n_genes_by_counts'], stat="count", bins=100, color='orange', kde=True, ax=ax2)
            plt.axvline(percentiles['n_genes_by_counts'][5], color='blue', linestyle='--')
            plt.axvline(percentiles['n_genes_by_counts'][95], color='blue', linestyle='--')
            ax2.set_xlim([0., 10000.])

            plt.tight_layout()
            plt.savefig(os.path.join(args.results, f'QC_Density_{sample}.png'))
            plt.close()
            
            print("\nVisualize density plot showing fraction of mitochondrial and ribosomal genes in {}".format(sample))
            ax1 = plt.subplot(1, 2, 1)
            sns.histplot(adata[adata.obs['sample']== sample].obs['pct_counts_mt'], stat="count", bins=100, kde=True, color='limegreen', ax=ax1)
            plt.axvline(MT_PERCENTAGE, 0, 1, c='red', linestyle='--')
            ax1.set_xlim([0., 25.])
            
            ax2 = plt.subplot(1, 2, 2)
            sns.histplot(adata[adata.obs['sample']== sample].obs['pct_counts_ribo'], stat="count", bins=100, kde=True, color='deepskyblue', ax=ax2)
            ax2.set_xlim([0., 60.])
            
            plt.savefig(os.path.join(args.results, f'QC_Density_MT-Ribo_{sample}.png'))
            plt.close()
# --------------------------------------------------------------------------------------------------------------------
#                           APPLY QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    # Filter cells of low quality

    print("\n===== FILTER CELLS BASED ON QUALITY METRICS =====")
    print('''Filter low quality cells on the basis of number of counts per barcode (count depth),number of genes per barcode of mitochondrial, and fraction of counts from mitochondrial genes per barcode    ''')

    #Filter based on MIN_COUNT
    sc.pp.filter_cells(adata, min_counts=percentiles['total_counts'][5],inplace=True)

    #Filter based on MIN_GENES
    sc.pp.filter_cells(adata, min_genes=percentiles['n_genes_by_counts'][5],inplace=True)

    print("Count matrix for combined samples has {} cells and {} genes after filtering".format(adata.shape[0],adata.shape[1]))
    #Filter based on MT_PERCENTAGE
    cell_number =adata[adata.obs.pct_counts_mt >= MT_PERCENTAGE].shape[0]
    print('''filter out {} cells for which the expression of mithocondrial genes is more than {}%'''.format(cell_number,MT_PERCENTAGE))
    adata = adata[adata.obs.pct_counts_mt < MT_PERCENTAGE, :]

    #Filter based on number of cells
    MIN_CELLS = round(1/100 * adata.shape[0])

    #Filter based on MIN_CELLS
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS, inplace=True)

    print("Count matrix has {} cells and {} genes".format(adata.shape[0],adata.shape[1]))


# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    
    fig, ax = plt.subplots(figsize=(20,10))
    print("\nVisualized the number of cells after filtering for each sample")
    sns.histplot(adata.obs, x="sample", stat="count", ax=ax)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90.)
    plt.savefig(os.path.join(args.results,'Cells_after_filtering.png'))
    plt.close()
    

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    
    print("\n===== SAVING OUTPUT FILE =====")

    print("Saving h5ad data to file {}".format(output))
    adata.write(output)
    print("Done!")
            
#####################################################################################################

if __name__ == '__main__':
    main()
