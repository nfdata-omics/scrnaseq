#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import argparse
import pathlib
import warnings
import os
import pandas as pd
import polars as pl
import numpy as np
import plotly.subplots as sp
import snapatac2 as snap
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image

warnings.filterwarnings("ignore")
# PARAMETERS
# set script version number
VERSION = "0.0.1"

import os
os.environ['SNAPATAC2_CACHE_DIR'] = '/home/camilla.callierotti/.snapatac2'

# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function calculates common quality control (QC) metrics for each sample and ATAC modality
    ,inspects QC plots for each computed QC metrics in each sample and filter cells based on QC metrics.
    """

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

#Define command line arguments with argparse
    parser = argparse.ArgumentParser(prog='QC_filter', usage='%(prog)s [options]', description = "QC metrics and filtering",
                                    epilog = "This function calculates common quality control (QC) metrics for each sample and ATAC modality, inspects QC plots for each sample and filters cells based on QC plots.",
                                    )
    parser.add_argument('-id', '--input-run-id', metavar='INPUT_RUN_ID', dest='input_run_id', nargs='+',
                        help="names of the run-id corresponding to the input adata")
    parser.add_argument('-fr','--input-fragment-files',metavar= 'FRAGMNET_FILES', type=pathlib.Path, nargs='+', dest='input_fragment_files',
                        required=True, help="paths of existing fragment file in tsv format")
    parser.add_argument('-fri','--input-fragment-files-index',metavar= 'FRAGMNET_FILES_INDEX', type=pathlib.Path,nargs='+', dest='input_fragment_files_index',
                        required=True, help="paths of existing index fragment file in tsv format")
    parser.add_argument('-t', '--tss_filter', dest='tss_threshold', type=float,default=2, help="parameters used to filter cells based on TSS score")
    parser.add_argument('-mif', '--min-fragments-counts', dest='min_fragments_counts', type=int, default=1000,
                        help="minimum number of fragments per cell to keep (default is 1000)")
    parser.add_argument('-maf', '--max-fragments-counts', dest='max_fragments_counts', type=int, default=100000,
                        help="maximum number of fragments per cell to keep (default is 100000)")
    parser.add_argument('-b', '--blacklist', metavar='BLACKLIST_FILE', type=pathlib.Path, default=None,
                        help="path to the blacklist file in bed format (default is None, no blacklist will be applied)")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.filtered_atac.h5ad",
                        help="path and name of the output h5ad file")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    input_run_id = args.input_run_id
    input_fragment_files = [str(f) for f in args.input_fragment_files]  
    input_fragment_files_index = [str(f) for f in args.input_fragment_files_index]
    output =args.out
    results_dir = args.results
    tss_threshold = args.tss_threshold
    min_fragments_counts = args.min_fragments_counts
    max_fragments_counts = args.max_fragments_counts
    blacklist_path = args.blacklist


    print("\n===== INPUT FRAGMENT FILES =====")
    for run, fragment in zip(input_run_id, input_fragment_files):
        print(f"Run: {run} - File: {fragment}")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ FRAGMENT FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the fragment files and compute basic QC metrics like the number of unique fragments per cell, fraction of duplicated reads and fraction of mitochondrial read
    print("\n===== READING FRAGMENT FILES =====")
    # Read the fragment file for each sample
    fragment_files = [str(f) for f in input_fragment_files]

    files = list(zip(input_run_id, fragment_files))
    print(files)

    print("\n===== VERIFYING FRAGMENT FILE MATCHING =====")
    for run_id, frag_file, frag_index in zip(input_run_id, input_fragment_files, input_fragment_files_index):
        print(f"Sample: {run_id}")
        print(f"  Fragment file: {frag_file}")
        print(f"  Index file: {frag_index}")
        # Verify files exist
        if not pathlib.Path(frag_file).exists():
            raise FileNotFoundError(f"Fragment file not found: {frag_file}")
        if not pathlib.Path(frag_index).exists():
            raise FileNotFoundError(f"Index file not found: {frag_index}")
        print("  Files verified.")

    # Import fragment files into AnnData objects and compute basic QC metrics like the number of unique fragments per cell, fraction of duplicated reads and fraction of mitochondrial read
    adatas_atac = snap.pp.import_fragments(
        fragment_files,
        file=[name + '.h5ad' for name in input_run_id],
        chrom_sizes=snap.genome.hg38,
        chrM=['chrM', 'M'],
        sorted_by_barcode=False,
        n_jobs =5
    )
    print("Fragment file import complete.")

    print("Print number of cells before filtering:")
    for run_id, adata in zip(input_run_id, adatas_atac):
        print(f"Sample {run_id}: {adata.n_obs} cells")

# --------------------------------------------------------------------------------------------------------------------
#                          COMPUTE AND VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    # Compute the fragment size distribution for each sample

    print("\n===== COMPUTE QUALITY METRICS {} =====")
    print(f"\n# Calculate the fragment size distribution for {input_fragment_files}")
    snap.metrics.frag_size_distr(adatas_atac,inplace = True,add_key='frag_size_distr',max_recorded_size=1000)

    print(f"\n# Calculate the tss enrichment for {input_fragment_files}")
    # Compute the TSS score for each sample
    snap.metrics.tsse(adatas_atac, snap.genome.hg38,inplace = True)
    for i, adata in enumerate(adatas_atac):
        fig2 = snap.pl.tsse(adata, show=False)
        fig2.show()
        fig2.write_image(f"TSS_score_sample_{i}.png", width=1000, height=600)

    print(f"\n# Compute the FRIP score for all samples")
    snap.metrics.frip(adatas_atac,regions= {"peaks_frac": snap.datasets.cre_HEA()}, normalized=True, count_as_insertion=False, inplace=True, n_jobs=5)
    print("\n===== FRIP SCORES PER SAMPLE =====")
    for run_id, adata in zip(input_run_id, adatas_atac):
        frip_score = adata.obs['peaks_frac'].mean()
        print(f"Sample {run_id}: FRIP score = {frip_score:.4f} (mean across cells)")
        print(f"  Min: {adata.obs['peaks_frac'].min():.4f}, Max: {adata.obs['peaks_frac'].max():.4f}")

    
    print("\n===== PLOT QUALITY METRICS =====")
    ### Fragment Size Distribution PDF ###
    
    print("Generating Fragment Size Distribution PDF...")
    with PdfPages(results_dir / "FragSizeDist_all_samples.pdf") as pdf:
        for run_id, adata in zip(input_run_id, adatas_atac):
            distr = adata.uns['frag_size_distr']
            print(f"Fragment size distribution for sample {run_id}:")
            print(distr) 
            if distr is None:
                print(f"No frag_size_distr found for {run_id}, skipping")
                continue

            # Convert distr to numeric arrays
            if isinstance(distr, dict):
                x_vals = np.array([int(k) for k in distr.keys()])
                y_vals = np.array([float(v) for v in distr.values()])
            else:
                y_vals = np.array(distr, dtype=float)
                x_vals = np.arange(len(y_vals))

            fig, ax = plt.subplots(figsize=(14, 8))

            ax.bar(x_vals, y_vals, color='chocolate', width=1)
            ax.set_xlim(left=0)
            ax.set_xlabel("Fragment Size")
            ax.set_ylabel("Count")
            ax.ticklabel_format(style='plain', axis='y')
            ax.grid(True, linestyle='--', alpha=0.5)
            sns.despine(ax=ax)

            fig.suptitle(f"QC – Fragment Size Distribution\n\nSample: {run_id}",
                        fontsize=18, fontweight='bold', y=1.02)

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    print("Fragment Size Distribution PDF generated.")
            
    from matplotlib.colors import LinearSegmentedColormap
    colors = ["#e1edf8", "#cbdff1", "#cbdff1", "#abd0e6",
            "#82badb", "#59a2cf", "#3787c0", "#1b6aaf",
            "#074d97", "#07306a"]
    blue_cmap = LinearSegmentedColormap.from_list("blue_kde", colors)

    print("Generating TSS Enrichment Score PDF...")
    with PdfPages(results_dir / "TSSE_score_all_samples.pdf") as pdf:
        for run_id, adata in zip(input_run_id, adatas_atac):
            png_file = f"tsse_{run_id}.png"
            # Increase width and height for higher resolution
            snap.pl.tsse(adata, out_file=png_file, show=False, width=1400, height=800)
            
            img = Image.open(png_file)
            
            fig, ax = plt.subplots(figsize=(14, 8))
            
            ax.imshow(img)
            ax.axis('off')
            
            fig.suptitle(f"QC – TSS Enrichment vs Unique Fragments\n\nSample: {run_id}",
                        fontsize=18, fontweight='bold', y=1.02)
            
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            
            print(f"Added TSSE plot for {run_id}")
            
    print("TSS Enrichment Score PDF generated.")
    
    print("\n===== PLOT QC HISTOGRAMS (CELL NUMBERS) =====")
    
    with PdfPages(results_dir / "QC_Histograms_all_samples.pdf") as pdf:
        for run_id, adata in zip(input_run_id, adatas_atac):

            n_fragment = adata.obs['n_fragment'].drop_nans()
            frac_dup = adata.obs['frac_dup'].drop_nans()
            peaks_frac = adata.obs['peaks_frac'].drop_nans()

            fig, axs = plt.subplots(1, 3, figsize=(21, 7))

            frag_counts, frag_bins = np.unique(n_fragment, return_counts=True)
            axs[0].bar(frag_counts, frag_bins, color='chocolate', width=1)
            axs[0].set_xscale("log")
            axs[0].set_xlabel("Number of Fragments")
            axs[0].set_ylabel("Cell Count")
            axs[0].ticklabel_format(style='plain', axis='y')
            axs[0].grid(True, linestyle='--', alpha=0.5)
            sns.despine(ax=axs[0])
            axs[0].set_title("Fragments per Barcode Distribution", fontsize=14, fontweight='bold')

            sns.histplot(frac_dup, bins=100, color="orange", ax=axs[1])
            axs[1].set_xlim(left=0)
            axs[1].set_xlabel("Fraction Duplicated")
            axs[1].set_ylabel("Cell Count")
            axs[1].grid(True, linestyle='--', alpha=0.5)
            sns.despine(ax=axs[1])
            axs[1].set_title("Fraction Duplicated per cell", fontsize=14, fontweight='bold')

            sns.histplot(peaks_frac, bins=100, color="deepskyblue", ax=axs[2])
            axs[2].set_xlabel("FRIP Score")
            axs[2].set_ylabel("Cell Count")
            axs[2].grid(True, linestyle='--', alpha=0.5)
            sns.despine(ax=axs[2])
            axs[2].set_title("FRIP Score per cell", fontsize=14, fontweight='bold')

            fig.suptitle(f"QC Metrics\n\nSample: {run_id}", fontsize=18, fontweight='bold', y=1.03)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    # Compute summary statistics by chromosome
    snap.metrics.summary_by_chrom(adatas_atac, mode='sum')
        
# --------------------------------------------------------------------------------------------------------------------
#                           FILTERS CELLS BASED ON QC METRICS
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== FILTER CELLS BASED ON QC METRICS =====")
    # Filter cells based on the TSS score and nucleosome signal
    snap.pp.filter_cells(adatas_atac, min_counts=min_fragments_counts, max_counts=max_fragments_counts,min_tsse=tss_threshold)

    print("Print number of cells after filtering for QC:")
    for run_id, adata in zip(input_run_id, adatas_atac):
        print(f"Sample {run_id}: {adata.n_obs} cells")

# --------------------------------------------------------------------------------------------------------------------
#                           CELL BY BIN MATRIX and DOUBLETS DETECTION
# --------------------------------------------------------------------------------------------------------------------

    # Save number of cells before filtering
    cell_counts = pd.DataFrame({"sample": input_run_id})
    cells_before_filtering = [adata.n_obs for adata in adatas_atac]
    cell_counts["cells_before_filtering"] = cells_before_filtering
    
    print("\n===== CELL BY BIN MATRIX =====")
    # Compute the tile matrix for each sample with a bin size of 500 bp
    snap.pp.add_tile_matrix(adatas_atac, bin_size=500, exclude_chroms=["chrM"], min_frag_size=None, max_frag_size=None, counting_strategy='paired-insertion', inplace=True)

    # Identify the most variable features in each sample based on the tile matrix
    print("\n===== IDENTIFY MOST VARIABLE FEATURES =====")
    snap.pp.select_features(adatas_atac, blacklist=blacklist_path, inplace=True)
    
    # Compute the doublet score for each sample
    print("\n===== COMPUTE DOUBLETS SCORE =====")
    snap.pp.scrublet(adatas_atac, n_comps=10, features = 'selected', sim_doublet_ratio=1,expected_doublet_rate = 0.10, n_jobs=2, use_approx_neighbors=True, inplace = True)

    # Filter out doublets based on the doublet score
    print("\n===== FILTER OUT DOUBLETS =====")
    snap.pp.filter_doublets(adatas_atac,verbose=True)

    # Save number of cells after doublet filtering
    cells_after_doublets = [adata.n_obs for adata in adatas_atac]
    cell_counts["cells_after_doublets"] = cells_after_doublets
    
    # Create csv of cell counts once for safety     
    cell_counts = cell_counts.sort_values("sample")
    cell_counts.to_csv(results_dir / "cell_counts_filters.csv",
                       index=False)

    print("Print number of cells after filtering doublets:")
    for run_id, adata in zip(input_run_id, adatas_atac):
        print(f"Sample {run_id}: {adata.n_obs} cells ") 


# --------------------------------------------------------------------------------------------------------------------
#                           CREATE ANNDATASET OBJECT
# --------------------------------------------------------------------------------------------------------------------

    # Create an AnnDataSet object to store the processed data
    # print("\n===== CREATE ANNDATASET OBJECT =====")
    data = snap.AnnDataSet(
    adatas=[(name, adata) for (name, _), adata in zip(files, adatas_atac)],
    filename=output
    )

    unique_cell_ids = [sa + ':' + bc for sa, bc in zip(data.obs['sample'], data.obs_names)]
    data.obs_names = unique_cell_ids
    assert data.n_obs == np.unique(data.obs_names).size

    print("\nNumber of cells per sample after filtering:")
    print(data.obs['sample'].value_counts())
    
    print("\n===== SAVE OBS INTO ANNDATASET =====")
    data.obs['tsse'] = data.adatas.obs['tsse']
    data.obs['n_fragment'] = data.adatas.obs['n_fragment']
    data.obs['frac_dup'] = data.adatas.obs['frac_dup']
    data.obs['frac_mito'] = data.adatas.obs['frac_mito']
    data.obs['peaks_frac'] = data.adatas.obs['peaks_frac']
    data.obs['doublet_score'] = data.adatas.obs['doublet_score']
    data.obs['doublet_probability'] = data.adatas.obs['doublet_probability']
    data.obsm['fragment_paired'] = data.adatas.obsm['fragment_paired']
    
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    data.close()

    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
