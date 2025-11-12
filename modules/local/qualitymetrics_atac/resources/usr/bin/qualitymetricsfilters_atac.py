#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import argparse
import pathlib
import warnings

import numpy as np
import plotly.subplots as sp
import snapatac2 as snap

warnings.filterwarnings("ignore")
# PARAMETERS
# set script version number
VERSION = "0.0.1"


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
    parser.add_argument('-n', '--nucleosome_filter', dest='nucleosome_threshold',type=float, default=2, help="parameters used to filter cells based on nucleosome signal")
    parser.add_argument('-t', '--tss_filter', dest='tss_threshold', type=float,default=1, help="parameters used to filter cells based on TSS score")
    parser.add_argument('-mif', '--min-fragments-counts', dest='min_fragments_counts', type=int, default=5000,
                        help="minimum number of fragments per cell to keep (default is 5000)")
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
    nucleosome_threshold = args.nucleosome_threshold
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

    # Import fragment files into AnnData objects and compute basic QC metrics like the number of unique fragments per cell, fraction of duplicated reads and fraction of mitochondrial read
    adatas_atac = snap.pp.import_fragments(
        fragment_files,
        file=[name + '.h5ad' for name in input_run_id],
        chrom_sizes=snap.genome.hg38,
        chrM=['chrM', 'M'],
        sorted_by_barcode=False,
        n_jobs =5
    )
    print(adatas_atac)

# --------------------------------------------------------------------------------------------------------------------
#                          COMPUTE AND VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    # Compute the fragment size distribution for each sample

    print("\n===== COMPUTE QUALITY METRICS {} =====")
    print(f"\n# Calculate the fragment size distribution for {input_fragment_files}")
    snap.metrics.frag_size_distr(adatas_atac,inplace = True,add_key='frag_size_distr',max_recorded_size=1000)
    
    for i, adata in enumerate(adatas_atac):
        fig1 = snap.pl.frag_size_distr(adata, show=False)
        fig1.show()
        fig1.update_yaxes(type="log")
        fig1.write_image(f"FragSizeDist_sample_{i}.png", width=1000, height=600)
        

    print(f"\n# Compute TSS enrichment score per cell for sample {input_fragment_files}")
    # Compute the TSS score for each sample
    snap.metrics.tsse(adatas_atac, snap.genome.hg38,inplace = True)
    for i, adata in enumerate(adatas_atac):
        fig2 = snap.pl.tsse(adata, show=False)
        fig2.show()
        fig2.write_image(f"TSS_score_sample_{i}.png", width=1000, height=600)

    print(f"\n# Compute the FRIP score for all samples")
    snap.metrics.frip(adatas_atac,regions= {"peaks_frac": snap.datasets.cre_HEA()}, normalized=True, count_as_insertion=False, inplace=True, n_jobs=5)

    # Compute summary statistics by chromosome
    snap.metrics.summary_by_chrom(adatas_atac, mode='sum')
    
    print(adatas_atac)
    
# --------------------------------------------------------------------------------------------------------------------
#                           FILTERS CELLS BASED ON QC METRICS
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== FILTER CELLS BASED ON QC METRICS =====")
    # Filter cells based on the TSS score and nucleosome signal
    snap.pp.filter_cells(adatas_atac, min_counts=min_fragments_counts, max_counts=max_fragments_counts,min_tsse=tss_threshold)
    #print(f"After filtering, MuData matrix for combined samples has {adatas_atac.shape[0]} cells and {adatas_atac.shape[1]} fragments")
# --
# --------------------------------------------------------------------------------------------------------------------
#                           CELL BY BIN MATRIX and DOUBLETS DETECTION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== CELL BY BIN MATRIX =====")
    # Compute the tile matrix for each sample with a bin size of 500 bp
    snap.pp.add_tile_matrix(adatas_atac, bin_size=500, exclude_chroms=["chrM"], min_frag_size=None, max_frag_size=None, counting_strategy='paired-insertion', inplace=True)

    # Identify the most variable features in each sample based on the tile matrix
    print("\n===== IDENTIFY MOST VARIABLE FEATURES =====")
    snap.pp.select_features(adatas_atac, blacklist=blacklist_path, inplace=True)
    
    # Compute the doublet score for each sample
    print("\n===== COMPUTE DOUBLETS SCORE =====")
    snap.pp.scrublet(adatas_atac,inplace = True)

    # Filter out doublets based on the doublet score
    print("\n===== FILTER OUT DOUBLETS =====")
    snap.pp.filter_doublets(adatas_atac)

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
