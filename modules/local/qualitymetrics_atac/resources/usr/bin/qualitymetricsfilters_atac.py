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
import snapatac2 as snap
import glob
import anndata as ad





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
    parser.add_argument('-fr','--input-fragment-files',metavar= 'FRAGMNET_FILES',type=pathlib.Path, nargs='+',dest='input_fragment_files',
                        required=True, help="paths of existing fragment file in tsv format")
    parser.add_argument('-fri','--input-fragment-files-index',metavar= 'FRAGMNET_FILES_INDEX',type=pathlib.Path,nargs='+', dest='input_fragment_files_index',
                        required=True, help="paths of existing index fragment file in tsv format")
    parser.add_argument('-n', '--nucleosome_filter',dest='nucleosome_threshold',type=float,default=2,help="parameters used to filter cells based on nucleosome signal")
    parser.add_argument('-t', '--tss_filter',dest='tss_threshold',type=float,default=1,help="parameters used to filter cells based on TSS score")
    parser.add_argument('-c', '--min-fragments-counts', dest='min_fragments_counts', type=int, default=1000,
                        help="minimum number of fragments per cell to keep (default is 1000)")
    parser.add_argument('-C', '--max-fragments-counts', dest='max_fragments_counts', type=int, default=100000,
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

    print("\n===== INPUT H5AD FILES =====")
    input_run_id = args.input_run_id
    input_fragment_file = str(args.input_fragment_files)
    input_fragment_file_index = str(args.input_fragment_files_index)
    output =args.out
    nucleosome_threshold = args.nucleosome_threshold
    tss_threshold = args.tss_threshold
    min_fragments_counts = args.min_fragments_counts
    max_fragments_counts = args.max_fragments_counts
    blacklist_path = args.blacklist


    for run, fragment in zip(input_run_id, input_fragment_file):
        print(f"Run: {run} - File: {fragment}")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ FRAGMENT FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the fragment files and compute basic QC metrics like the number of unique fragments per cell, fraction of duplicated reads and fraction of mitochondrial read
    print("\n===== READING FRAGMENT FILES =====")
    # read the fragment file for each sample
    print(f"\nProcessing fragment file object {input_fragment_file} ... ", end ='')

    fragment_files = []
    for folder in glob.glob("*/atac_fragments.tsv.gz"):
        fragment_files.append(folder)

    print(fragment_files)


    files = list(zip(input_run_id, fragment_files))
    print(files)

    adatas_atac = snap.pp.import_fragments(
        [fl for _, fl in files],
        file=[name + '.h5ad' for name, _ in files],
        chrom_sizes=snap.genome.hg38,
        chrM=['chrM', 'M'],
        sorted_by_barcode=False
    )

# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE AND VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------

    # Compute the fragment size distribution for each sample

    print("\n===== COMPUTE QUALITY METRICS {} =====")
    print(f"\n# Calculate the fragment size distribution for {input_fragment_file}")
    snap.metrics.frag_size_distr(adatas_atac,inplace = True,add_key='frag_size_distr',max_recorded_size=1000)
    for i, adata in enumerate(adatas_atac):
    #for adata in adatas_atac:
        #sample_name = adata.obs["sample"].unique()[0]
        fig1 = snap.pl.frag_size_distr(adata, show=False)  # genera la figura per il singolo campione
        fig1.show()
        fig1.update_yaxes(type="log")  # scala logaritmica sull'asse y
        fig1.write_image(f"FragSizeDist_sample_{i}.png", width=1000, height=600)
        #fig1.write_image(f"FragSizeDist_{sample_name}.png", width=1000, height=600) #To save with sample name


    print(f"\n# Calculate the nucleosome signal for {input_fragment_file}")
    # Compute the nucleosome signal for each sample
    '''
    frag_dist = adatas_atac.uns['frag_size_distr']
    nuc_signal = []

    for cell_id in adatas_atac.obs_names:
        dist = frag_dist[cell_id]  # assicurati che frag_dist[cell_id] dia un array di counts
        subnuc = dist[38:91].sum()
        mono = dist[180:248].sum()
        if subnuc == 0:
            nuc_signal.append(np.nan)
        else:
            nuc_signal.append(mono / subnuc)

    adatas_atac.obs['nucleosome_signal'] = nuc_signal
    '''

    print(f"\n# Compute TSS enrichment score per cell for sample {input_fragment_file}")
    # Compute the TSS score for each sample
    snap.metrics.tsse(adatas_atac, snap.genome.hg38,inplace = True)
    for i, adata in enumerate(adatas_atac):
        fig2 = snap.pl.tsse(adata, show=False)
        fig2.show()
        fig2.write_image(f"TSS_score_sample_{i}.png", width=1000, height=600)


    print(f"\n# Calculate the FRIP for {input_fragment_file}")
    snap.metrics.frip(adatas_atac,{"peaks_frac": snap.datasets.cre_HEA()},normalized=True,inplace=True)

    print(f"\n# Calculate the metric summary for each chrom for {input_fragment_file}")
    snap.metrics.summary_by_chrom(adatas_atac, mode='sum')



# --------------------------------------------------------------------------------------------------------------------
#                           FILTERS CELLS BASED ON QC METRICS
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== FILTER CELLS BASED ON QC METRICS =====")
    # Filter cells based on the TSS score and nucleosome signal
    snap.pp.filter_cells(adatas_atac, min_counts=min_fragments_counts, max_counts=max_fragments_counts,min_tsse=tss_threshold)
    #filtered_adatas = []
    #for adata in adatas_atac:
    #    adata_filtered = adata[(adata.obs['frac_dup'] <= 0.2) & (adata.obs['frac_mito'] <= 0.3)].copy()
    #    filtered_adatas.append(adata_filtered)

    #adatas_atac = filtered_adatas
# --------------------------------------------------------------------------------------------------------------------
#                           CELL BY BIN MATRIX and DOUBLETS DETECTION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== CELL BY BIN MATRIX =====")
    # Compute the tile matrix for each sample with a bin size of 500 bp
    snap.pp.add_tile_matrix(adatas_atac,bin_size=500,exclude_chroms=["chrM"],min_frag_size=50,max_frag_size=1000,counting_strategy='paired-insertion',inplace=True)

    # Identify the most variable features in each sample based on the tile matrix
    print("\n===== IDENTIFY MOST VARIABLE FEATURES =====")
    #snap.pp.select_features(adatas_atac, n_features=500000,inplace = True)
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
    print("\n===== CREATE ANNDATASET OBJECT =====")
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

    '''
    frag_sizes = []
    nucleosome_signals = []
    library_tsse = []
    frac_overlap_TSS = []


    for name, adata in data.adatas:
        if 'frag_size_distr' in adata.uns:
            frag_sizes.append(pd.Series(
                adata.uns['frag_size_distr'], index=adata.obs.index
            ))
        if 'nucleosome_signal' in adata.uns:
            nucleosome_signals.append(pd.Series(
                adata.uns['nucleosome_signal'], index=adata.obs.index
            ))
        if 'library_tsse' in adata.uns:
            library_tsse.append(pd.Series(
                adata.uns['library_tsse'], index=adata.obs.index
            ))
        if 'frac_overlap_TSS' in adata.uns:
            frac_overlap_TSS.append(pd.Series(
                adata.uns['frac_overlap_TSS'], index=adata.obs.index
            ))
    if frag_sizes:
        data.obs['frag_size_distr'] = pd.concat(frag_sizes)
    if nucleosome_signals:
        data.obs['nucleosome_signal'] = pd.concat(nucleosome_signals)
    if library_tsse:
        data.obs['library_tsse'] = pd.concat(library_tsse)
    if frac_overlap_TSS:
        data.obs['frac_overlap_TSS'] = pd.concat(frac_overlap_TSS)
    '''
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    data.close()


    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
