#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import argparse                     # command line arguments parser
import warnings
import os                           # filesystem utilities
import pathlib
from pathlib import Path                      # library for handle filesystem paths
import numpy as np
import pandas as pd                 # library for data analysis and manipulation
import scanpy as sc                 # single-cell data processing
import matplotlib.pyplot as plt     # library for visualization
import seaborn as sns               # library for statistical data visualization
import mudata as md
import muon as mu
from muon import atac as ac
from scipy.stats import median_abs_deviation
from matplotlib.backends.backend_pdf import PdfPages


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
                                    epilog = "This function calculates common quality control (QC) metrics for each sample and modality, inspects QC plots for each sample and filters cells based on QC plots.")
    parser.add_argument('-ad','--input-h5mu-combined', metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing matrix files in h5mu format (including file names)")
    parser.add_argument('-d','--input-csv-doublets', metavar= 'CSV_DOUBLETS_TABLE', type=pathlib.Path, dest='input_csv_table',
                        default=pathlib.Path(''), help="paths of existing doublets table in csv format")
    parser.add_argument('-mt', '--mt-thresold', dest='mt_threshold', type=float, default=15, help="parameters used to filter cells based on mithocondrial gene content")
    parser.add_argument('-min', '--min-umi', dest='min_umi_gex', type=int, default=1400,
                        help="minimum number of UMI per cell to keep (default is 1400)")
    parser.add_argument('-max', '--max-umi', dest='max_umi_gex', type=int, default=15000,
                        help="maximum number of UMI per cell to keep (default is 15000)")
    parser.add_argument('-ming', '--min-genes', dest='min_genes_gex', type=int, default=200,
                        help="minimum number of genes per cell to keep (default is 200)")
    parser.add_argument('-maxg', '--max-genes', dest='max_genes_gex', type=int, default=5000,
                        help="maximum number of genes per cell to keep (default is 5000)")
    parser.add_argument('-minc', '--min-cells', dest='min_cells_gex', type=int, default=5,
                        help="minimum number of cells per gene to keep (default is 5)")
    parser.add_argument('-minf', '--min-features-adt', dest='min_features_adt', type=int, default=3,
                        help="minimum number of features per cell to keep (default is 3)")
    parser.add_argument('-mincadt', '--min-counts-adt', dest='min_counts_adt', type=int, default=500,
                        help="minimum number of counts per cell to keep (default is 500)")
    parser.add_argument('-csv', '--csv_out', metavar='QUALITY_CONTROL', default="summary_qualitycontrol.csv",
                        help="path and name of csv table with ranked marker genes for each cluster and resolution")
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
    output_csv= Path(args.csv_out)
    output =args.out
    mt_threshold = args.mt_threshold
    min_umi_gex = args.min_umi_gex
    max_umi_gex = args.max_umi_gex
    min_genes_gex = args.min_genes_gex
    max_genes_gex = args.max_genes_gex
    min_cells_gex = args.min_cells_gex
    min_features_adt = args.min_features_adt
    min_counts_adt = args.min_counts_adt


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
        gex.var["gene_name_upper"] = gex.var["gene_name"].str.upper()

        gex.var["mt"] = gex.var["gene_name_upper"].str.startswith("MT-")
        gex.var["ribo"] = gex.var["gene_name_upper"].str.startswith(("RPS", "RPL"))
        gex.var["hb"] = gex.var["gene_name_upper"].str.startswith("HB") & ~gex.var["gene_name_upper"].str.startswith("HBP")
        sc.pp.calculate_qc_metrics(gex, qc_vars=["mt", "ribo", "hb"],percent_top=[20], log1p=True, inplace=True)
        print(gex.obs[['pct_counts_mt', 'pct_counts_ribo']].head())


# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------
# Visualize quality metrics: highest expressed genes, number of genes expressed, total counts per cell and fraction of mitochondrial, ribosomal and hemoglobin genes

        fig, ax = plt.subplots(figsize=(35,25))
        print("\nVisualized the number of cells for each pool before filtering")
        sns.histplot(gex.obs, x="sample", stat="count", ax=ax)
        locs, labels = plt.xticks()
        ax.set_xlabel("Sample name", fontsize=30)
        ax.set_ylabel("Cell number", fontsize=30)
        plt.setp(labels, rotation=90.,fontsize=30)
        ax.tick_params(axis='y', labelsize=30)
        plt.savefig(os.path.join(args.results,'Cells_before_filtering.pdf'), bbox_inches='tight', dpi=300)
        plt.close()

        with PdfPages(os.path.join(args.results, "QC_Density_all_samples.pdf")) as pdf_counts, PdfPages(os.path.join(args.results, "QC_Density_MT-Ribo_all_samples.pdf")) as pdf_mt_ribo:
            for sample in gex.obs['sample'].unique():

                # Counts distribution
                print(f"\nVisualizing density plot showing number of genes expressed, total counts per cell in {sample}")
                fig, axs = plt.subplots(1, 2, figsize=(14, 8))
                sns.histplot(gex[gex.obs['sample'] == sample].obs['total_counts'], stat="count", bins=500, color='chocolate', kde=True, ax=axs[0])
                axs[0].axvline(min_umi_gex, color='blue', linestyle='--')
                axs[0].axvline(max_umi_gex, color='blue', linestyle='--')
                axs[0].set_title("Total Counts per Cell")

                sns.histplot(gex[gex.obs['sample'] == sample].obs['n_genes_by_counts'], stat="count", bins=100, color='orange', kde=True, ax=axs[1])
                axs[1].axvline(min_genes_gex, color='blue', linestyle='--')
                axs[1].axvline(max_genes_gex, color='blue', linestyle='--')
                axs[1].set_title("Genes per Cell")

                fig.suptitle(f"QC – Counts and Genes per Cell\n\nSample: {sample}", fontsize=18, fontweight='bold', y=1.02)
                plt.tight_layout()
                pdf_counts.savefig(fig, bbox_inches='tight')
                plt.close(fig)

                # Mito/Ribo distribution
                print(f"\nVisualizing density plot showing fraction of mitochondrial and ribosomal genes in {sample}")
                fig, axs = plt.subplots(1, 2, figsize=(14, 8))
                sns.histplot(gex[gex.obs['sample'] == sample].obs['pct_counts_mt'], stat="count", bins=100, kde=True, color='limegreen', ax=axs[0])
                axs[0].axvline(mt_threshold, 0, 1, c='red', linestyle='--')
                axs[0].set_title("Mitochondrial Fraction")

                sns.histplot(gex[gex.obs['sample'] == sample].obs['pct_counts_ribo'], stat="count", bins=100, kde=True, color='deepskyblue', ax=axs[1])
                axs[1].set_title("Ribosomal Fraction")

                fig.suptitle(f"QC – Mitochondrial & Ribosomal Fractions\n\nSample: {sample}", fontsize=18, fontweight='bold', y=1.02)
                plt.tight_layout()
                pdf_mt_ribo.savefig(fig, bbox_inches='tight')
                plt.close(fig)

        print(f"\nSaved multi-page PDFs")


# --------------------------------------------------------------------------------------------------------------------
#                           EVALUATE GENES BASED ON NUMBER OF CELLS
# --------------------------------------------------------------------------------------------------------------------
        min_cells = round(1/100 * gex.shape[0])
        min_cells = max(min_cells, min_cells_gex)  # Ensure at least 3 cells for a gene to be considered expressed
        gex.var["gene_pass"] = gex.var["n_cells_by_counts"] >= min_cells

# --------------------------------------------------------------------------------------------------------------------
#                           EVALUATE CELLS BASED ON HARD FILTERS
# --------------------------------------------------------------------------------------------------------------------

        #Hard filtering based on total counts, number of genes expressed and fraction of mitochondrial genes
        gex.obs["total_counts_outlier"] = gex.obs["total_counts"] > max_umi_gex
        gex.obs["n_genes_by_counts_outlier"] = ((gex.obs["n_genes_by_counts"] < min_genes_gex) | (gex.obs["n_genes_by_counts"] > max_genes_gex))
        gex.obs["mt_outlier"] = is_outlier(gex, "pct_counts_mt", 3) | ( gex.obs["pct_counts_mt"] > mt_threshold )

        gex.obs["hard_filter_gex"] = (gex.obs["total_counts_outlier"] | gex.obs["n_genes_by_counts_outlier"] | gex.obs["mt_outlier"])

# --------------------------------------------------------------------------------------------------------------------
#                           EVALUATE CELLS BASED ON SOFT FILTERS
# --------------------------------------------------------------------------------------------------------------------
        print(gex.obs.columns)
        gex.obs["soft_filters_gex"] = (
        is_outlier(gex, "log1p_total_counts", 5)
        | is_outlier(gex, "log1p_n_genes_by_counts", 5)
        | is_outlier(gex, "pct_counts_in_top_20_genes", 5)
        | gex.obs["mt_outlier"]
        )

        hard_filter_gex_pool = gex.obs.groupby(['sample', 'hard_filter_gex']).size().unstack(fill_value=0).rename(columns={False: 'hard_filters_gex_pass', True: 'hard_filters_gex_fail'})
        soft_filter_gex_pool = gex.obs.groupby(['sample', 'soft_filters_gex']).size().unstack(fill_value=0).rename(columns={False: 'soft_filters_gex_pass', True: 'soft_filters_gex_fail'})
        #hard_filter_gex_sample = gex.obs.groupby(['Inferred_donor', 'hard_filter_gex',]).size().unstack(fill_value=0).rename(columns={False: 'hard_filters_gex_pass', True: 'hard_filters_gex_fail'})
        #soft_filter_gex_sample = gex.obs.groupby(['Inferred_donor', 'soft_filters_gex']).size().unstack(fill_value=0).rename(columns={False: 'soft_filters_gex_pass', True: 'soft_filters_gex_fail'})

# --------------------------------------------------------------------------------------------------------------------
#                           APPLY QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------


        print("\n===== FILTER CELLS BASED ON QUALITY METRICS =====")
        print('Filter low quality cells on the basis of number of counts per barcode (count depth),number of genes per barcode of mitochondrial, and fraction of counts from mitochondrial genes per barcode')

        #Filter based on MIN_COUNT
        mu.pp.filter_obs(gex, 'total_counts',lambda x: x >= min_umi_gex)

        #Filter based on MIN_GENES
        mu.pp.filter_obs(gex, 'n_genes_by_counts',lambda x: x >= min_genes_gex)

        print(f"Count matrix for combined samples has {gex.shape[0]} cells and {gex.shape[1]} genes after filtering")
        #Filter based on MT_PERCENTAGE
        cell_number =gex[gex.obs.pct_counts_mt >= mt_threshold].shape[0]
        print(cell_number)
        print(f'filter out {cell_number} cells for which the expression of mithocondrial genes is more than {mt_threshold}%')
        mu.pp.filter_obs(gex,'pct_counts_mt', lambda x: x < mt_threshold)

        #Filter based on number of cells
        mu.pp.filter_var(gex, 'n_cells_by_counts', lambda x: x >= min_cells)


        print(f"Count matrix has {gex.shape[0]} cells and {gex.shape[1]} genes")

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE FILTERING RESULT
# --------------------------------------------------------------------------------------------------------------------

        fig, ax = plt.subplots(figsize=(30,25))
        print("\nVisualized the number of cells after filtering for each sample")
        sns.histplot(gex.obs, x="sample", stat="count", ax=ax)
        locs, labels = plt.xticks()
        ax.set_xlabel("Sample name", fontsize=30)
        ax.set_ylabel("Cell number", fontsize=30)
        plt.setp(labels, rotation=90.,fontsize=30)
        ax.tick_params(axis='y', labelsize=30)
        plt.savefig(os.path.join(args.results,'Cells_after_filtering.pdf'), bbox_inches='tight', dpi=300)
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

        sc.pp.calculate_qc_metrics(pro, inplace=True,log1p=True,percent_top=[20])


# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------
# Visualize quality metrics: distribution of ADTs
# Added filters on samples with Antibody Capture feature type

        fig, ax = plt.subplots(figsize=(50,20))
        print("\nVisualized the distribution of ADTs per cell over all samples before filtering")
        for sample in pro.obs['sample'].unique():
            if (pro[pro.obs['sample'] == sample].var["feature_types"].str.contains("Antibody Capture")).any():
            #if (pro[pro.obs['sample'] == sample].var["feature_types"].str.contains("ab|ADT")).any():
                print(f"\nVisualized the distribution of ADTs per cell per {sample} before filtering ")
                ax1 = plt.subplot(1, 2, 1)
                sns.histplot(pro[pro.obs['sample']== sample].obs.total_counts,ax=ax1)
                plt.axvline(min_counts_adt, color='blue', linestyle='--')
                ax1.set_xlim([0., 8000.])
                ax2 = plt.subplot(1, 2, 2)
                sns.histplot(pro[pro.obs['sample']== sample].obs.n_genes_by_counts,ax=ax2)
                plt.axvline(min_features_adt, color='blue', linestyle='--')
                ax2.set_xlim([0., 200.])
                plt.savefig(os.path.join(args.results, f'ADTs_Distribution_{sample}.png'))
                plt.close()
            else:
                print(f"\nNo Antibody Capture feature type in {sample}")

# --------------------------------------------------------------------------------------------------------------------
#                           EVALUATE CELLS BASED ON HARD FILTERS
# --------------------------------------------------------------------------------------------------------------------


        #Hard filtering based on total counts, number of genes expressed and fraction of mitochondrial genes
        pro.obs["total_counts_outlier"] = pro.obs["total_counts"] < min_counts_adt
        pro.obs["n_genes_by_counts_outlier"] = pro.obs["n_genes_by_counts"] < min_features_adt

        pro.obs["hard_filter_pro"] = (pro.obs["total_counts_outlier"] | pro.obs["n_genes_by_counts_outlier"])

# --------------------------------------------------------------------------------------------------------------------
#                           EVALUATE CELLS BASED ON SOFT FILTERS
# --------------------------------------------------------------------------------------------------------------------

        pro.obs["soft_filter_pro"] = (
        is_outlier(pro, "log1p_total_counts", 5)
        | is_outlier(pro, "log1p_n_genes_by_counts", 5)
        | is_outlier(pro, "pct_counts_in_top_20_genes", 5)
        )



        hard_filter_pro_pool = pro.obs.groupby(['sample', 'hard_filter_pro']).size().unstack(fill_value=0).rename(columns={False: 'hard_filters_pro_pass', True: 'hard_filters_pro_fail'})
        soft_filter_pro_pool = pro.obs.groupby(['sample', 'soft_filter_pro']).size().unstack(fill_value=0).rename(columns={False: 'soft_filters_pro_pass', True: 'soft_filters_pro_fail'})
        #hard_filter_pro_sample = pro.obs.groupby(['Inferred_donor', 'hard_filter_pro']).size().unstack(fill_value=0).rename(columns={False: 'hard_filters_pro_pass', True: 'hard_filters_pro_fail'})
        #soft_filter_pro_sample = pro.obs.groupby(['Inferred_donor', 'soft_filter_pro']).size().unstack(fill_value=0).rename(columns={False: 'soft_filters_pro_pass', True: 'soft_filters_pro_fail'})


# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
        print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
        mdata.mod['pro'] = pro
        mdata.update()

    else:
        print("CITE modality does not exist in mdata.mod.")

# --------------------------------------------------------------------------------------------------------------------
#                           ADDED QUALITY METRICS INTO ADATA.OBS AND PRINT SUMMARY TABLE
# --------------------------------------------------------------------------------------------------------------------

    output_csv_pool = output_csv.with_name(output_csv.stem + "_by_pool.csv")
    print("\n===== ADDED QUALITY METRICS INTO ADATA.OBS AND PRINT SUMMARY TABLE =====")
    dfs = []
    for name in ['hard_filter_gex_pool', 'soft_filter_gex_pool','hard_filter_pro_pool','soft_filter_pro_pool']:
        df = locals().get(name)
        if df is not None:
            dfs.append(df)
    if dfs:
        summary_table = pd.concat(dfs, axis=1).fillna(0).astype(int)
        summary_table.to_csv(output_csv_pool)
        print("Done!")

    else:
        print("No dataframe available for concatenation.")

    output_csv_sample = output_csv.with_name(output_csv.stem + "_by_sample.csv")
    print("\n===== ADDED QUALITY METRICS INTO ADATA.OBS AND PRINT SUMMARY TABLE =====")
    dfs = []
    for name in ['hard_filter_gex_sample', 'soft_filter_gex_sample','hard_filter_pro_sample','soft_filter_pro_sample']:
        df = locals().get(name)
        if df is not None:
            dfs.append(df)
    if dfs:
        summary_table = pd.concat(dfs, axis=1).fillna(0).astype(int)
        summary_table.to_csv(output_csv_sample)
        print("Done!")

    else:
        print("No dataframe available for concatenation.")


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
