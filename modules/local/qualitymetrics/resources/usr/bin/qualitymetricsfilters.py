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
import muon as mu
from muon import atac as ac
from scipy.stats import median_abs_deviation



warnings.filterwarnings("ignore")
# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAD FUNCTION
# ====================================================================================================================

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

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
                                    epilog = "This function calculates common quality control (QC) metrics for each sample and modality, inspects QC plots for each sample and filters cells based on QC plots.",
                                    )
    parser.add_argument('-ad','--input-h5mu-combined',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing matrix files in h5mu format (including file names)")
    parser.add_argument('-d','--input-csv-doublets',metavar= 'CSV_DOUBLETS_TABLE', type=pathlib.Path, dest='input_csv_table',
                        default=pathlib.Path(''),help="paths of existing doublets table in csv format")
    parser.add_argument('-mt', '--mt-thresold',dest='mt_threshold',type=float,default=15,help="parameters used to filter cells based on mithocondrial gene content")
    parser.add_argument('-csv', '--csv_out', metavar='QUALITY_CONTROL', default="summary_qualitycontrol.csv",
                        help="path and name of excel table with ranked marker genes for each cluster and resolution")
    parser.add_argument('-o', '--out', metavar='H5MU_OUTPUT_FILE', type=pathlib.Path, default="matrix.filtered.h5mu",
                        help="path and name of the output h5mu file")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5MU FILES =====")
    input_h5mu_file = args.input_h5mu_files
    input_csv_table = args.input_csv_table
    output_csv= args.csv_out
    output =args.out
    mt_threshold = args.mt_threshold


    # print info on the available matrices
    print("Reading combined matrix from the following file:")
    print(f"-File {input_h5mu_file}")
    
# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5MU FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print(f"\nProcessing MuData object in folder {input_h5mu_file} ... ", end ='')

    mdata= md.read(input_h5mu_file)
    print("Done!")
    print(f"MuData matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== GEX MODALITY DATA =====")
    # Check if 'gex' exists in mdata.mod
    if 'gex' in mdata.mod:
        gex = mdata.mod['gex']

# --------------------------------------------------------------------------------------------------------------------
#                                 FILTER DOUBLETS
# --------------------------------------------------------------------------------------------------------------------
        #print("\n===== READING DOUBLETS TABLE =====")
        if input_csv_table and input_csv_table != pathlib.Path(''):
            input_csv_table=pd.read_csv(input_csv_table,index_col=0)

            gex.obs["doublets"] = input_csv_table['scDblFinder.class']
            gex.obs["doublet_score"] = input_csv_table['scDblFinder.score']
            #print("\n===== FILTER OUT DOUBLETS =====")
            #Filter based on doublet score
            cell_doublets =gex[gex.obs["doublets"] == 'doublet'].shape[0]
            print(f'''filter out {cell_doublets} cells which are doublets''')
            gex = gex[gex.obs.doublets != 'doublet', :]
            del gex.obs["doublets"]
            del gex.obs["doublet_score"]

# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

        # Compute the fraction of mitochondrial, ribosomal and hemoglobin genes

        print("\n===== COMPUTE QUALITY METRICS {} =====")
        print(f"\nCompute fraction of mitochondrial, ribosomal and hemoglobin genes for {input_h5mu_file}")
        gex.var["mt"] = gex.var["gene_symbols"].str.startswith("MT-")
        gex.var["ribo"] = gex.var["gene_symbols"].str.startswith(("RPS", "RPL")) 
        gex.var["hb"] = gex.var["gene_symbols"].str.startswith(("^HB[^(P)]"))
        sc.pp.calculate_qc_metrics(gex, qc_vars=["mt", "ribo", "hb"],percent_top=[20], log1p=True, inplace=True)
        print(gex.obs[['pct_counts_mt', 'pct_counts_ribo']].head())

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
#                           EVALUATE CELLS BASED ON MAD
# --------------------------------------------------------------------------------------------------------------------
        print(gex.obs.columns)
        gex.obs["MAD_outlier"] = (
        is_outlier(gex, "log1p_total_counts", 5)
        | is_outlier(gex, "log1p_n_genes_by_counts", 5)
        | is_outlier(gex, "pct_counts_in_top_20_genes", 5)
        )
        

# --------------------------------------------------------------------------------------------------------------------
#                           ADDED QUALITY METRICS INTO ADATA.OBS AND PRINT SUMMARY TABLE
# --------------------------------------------------------------------------------------------------------------------

        gex.obs["total_counts_outlier"] = ((gex.obs["total_counts"] < percentiles['total_counts'][5]) | (gex.obs["total_counts"] > percentiles['total_counts'][95]))
        gex.obs["n_genes_by_counts_outlier"] = ((gex.obs["n_genes_by_counts"] < percentiles['n_genes_by_counts'][5]) | (gex.obs["n_genes_by_counts"] > percentiles['n_genes_by_counts'][95]))
        gex.obs["mt_outlier"] = (gex.obs["pct_counts_mt"] > mt_threshold)

        counts_df = gex.obs.groupby(['sample', 'total_counts_outlier']).size().unstack(fill_value=0).rename(columns={False: 'total_counts_pass', True: 'total_counts_fail'})
        genes_df = gex.obs.groupby(['sample', 'n_genes_by_counts_outlier']).size().unstack(fill_value=0).rename(columns={False: 'n_genes_pass', True: 'n_genes_fail'})
        mt_df     = gex.obs.groupby(['sample', 'mt_outlier']).size().unstack(fill_value=0).rename(columns={False: 'mt_pass', True: 'mt_fail'})
        mad_df    = gex.obs.groupby(['sample', 'MAD_outlier']).size().unstack(fill_value=0).rename(columns={False: 'MAD_pass', True: 'MAD_fail'})

        summary_table = pd.concat([counts_df, genes_df, mt_df,mad_df], axis=1)
        summary_table = summary_table.fillna(0).astype(int)
        summary_table.to_csv(output_csv)
        print("Done!")

        '''
        print("\nVisualized the number of outliers for each sample")
        for sample in gex.obs['sample'].unique()
            fig, ax = plt.subplots(figsize=(20,10))
            outlier_counts.plot(kind='bar', stacked=True, ax=ax)
            plt.title("Outliers based on total counts")
            plt.xlabel("Sample")
            plt.ylabel("Number of cells")
            plt.legend(title="Outlier", labels=["Not Outlier", "Outlier"])
            plt.savefig(os.path.join(args.results,'Outliers_total_counts.png'))
            plt.close()
            fig, ax = plt.subplots(figsize=(20,10))
            outlier_genes.plot(kind='bar', stacked=True, ax=ax)
            plt.title("Outliers based on number of genes")
            plt.xlabel("Sample")
            plt.ylabel("Number of cells")
            plt.legend(title="Outlier", labels=["Not Outlier", "Outlier"])
            plt.savefig(os.path.join(args.results,'Outliers_genes.png'))
            plt.close()
            fig, ax = plt.subplots(figsize=(20,10))
            outlier_mt.plot(kind='bar', stacked=True, ax=ax)
            plt.title("Outliers based on mitochondrial genes")
            plt.xlabel("Sample")
            plt.ylabel("Number of cells")
            plt.legend(title="Outlier", labels=["Not Outlier", "Outlier"])
            plt.savefig(os.path.join(args.results,'Outliers_mitochondrial.png'))
            plt.close()
        '''

# --------------------------------------------------------------------------------------------------------------------
#                           APPLY QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

        if input_csv_table and input_csv_table != pathlib.Path(''):
        # Filter cells of low quality

            print("\n===== FILTER CELLS BASED ON QUALITY METRICS =====")
            print('Filter low quality cells on the basis of number of counts per barcode (count depth),number of genes per barcode of mitochondrial, and fraction of counts from mitochondrial genes per barcode')

            #Filter based on MIN_COUNT
            mu.pp.filter_obs(gex, 'total_counts',lambda x: x >= percentiles['total_counts'][5])

            #Filter based on MIN_GENES
            mu.pp.filter_obs(gex, 'n_genes_by_counts',lambda x: x >= percentiles['n_genes_by_counts'][5])

            print(f"Count matrix for combined samples has {gex.shape[0]} cells and {gex.shape[1]} genes after filtering")
            #Filter based on MT_PERCENTAGE
            cell_number =gex[gex.obs.pct_counts_mt >= mt_threshold].shape[0]
            print(cell_number)
            print(f'filter out {cell_number} cells for which the expression of mithocondrial genes is more than {mt_threshold}%')
            mu.pp.filter_obs(gex,'pct_counts_mt', lambda x: x < mt_threshold)

            #Filter based on number of cells
            min_cells = round(1/100 * gex.shape[0])
            mu.pp.filter_var(gex, 'n_cells_by_counts', lambda x: x >= min_cells)


            print(f"Count matrix has {gex.shape[0]} cells and {gex.shape[1]} genes")

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
#                                 CITE MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== CITE MODALITY DATA =====")
    # Check if 'pro' exists in mdata.mod
    if 'pro' in mdata.mod:
        pro = mdata.mod['pro']


# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    # Evaluate the distribution of ADTs per cells over all samples....

        print("\n===== COMPUTE QUALITY METRICS {} =====")
        print(f"\nCompute the distribution of ADTs per cells over all samples in {input_h5mu_file}")

        sc.pp.calculate_qc_metrics(pro, inplace=True,log1p=False,percent_top=None)
        
# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------
# Visualize quality metrics: distribution of ADTs
# Added filters on samples with Antibody Capture feature type

        fig, ax = plt.subplots(figsize=(40,10))
        print("\nVisualized the distribution of ADTs per cell over all samples before filtering")
        for sample in pro.obs['sample'].unique():
            if (pro[pro.obs['sample'] == sample].var["feature_types"].str.contains("Antibody Capture")).any():
            #if (pro[pro.obs['sample'] == sample].var["feature_types"].str.contains("ab|ADT")).any():
                print(f"\nVisualized the distribution of ADTs per cell per {sample} before filtering ")
                ax1 = plt.subplot(1, 2, 1)
                sns.histplot(pro[pro.obs['sample']== sample].obs.total_counts,ax=ax1)
                ax1.set_xlim([0., 1000.])
                ax2 = plt.subplot(1, 2, 2)
                sns.histplot(pro[pro.obs['sample']== sample].obs.n_genes_by_counts,ax=ax2)
                ax2.set_xlim([0., 100.])
                plt.savefig(os.path.join(args.results, f'ADTs_Distribution_{sample}.png'))
                plt.close()
            else:
                print(f"\nNo Antibody Capture feature type in {sample}")
    else:
        print("CITE modality does not exist in mdata.mod.")

#
    
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
