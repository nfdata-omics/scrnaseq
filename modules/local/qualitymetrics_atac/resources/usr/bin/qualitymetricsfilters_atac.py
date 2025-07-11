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
    blacklist_path = args.blacklist


    # print info on the available matrices
    #print("Reading fragment file from the following file:")
    #print(f"-File {input_fragment_file}")

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
    #fig = snap.pl.frag_size_distr(adatas_atac, show=False)
    #fig.update_yaxes(type="log")
    #fig.write_image("Nucleosome_signal.png", width=1000, height=600)

    print(f"\n# Calculate the TSS score for {input_fragment_file}")
    # Compute the TSS score for each sample
    snap.metrics.tsse(adatas_atac, snap.genome.hg38,inplace = True)
    #fig2 = snap.pl.tsse(adatas_atac, interactive=False)
    #fig2.write_image("TSS_score.png", width=1000, height=600)

    print(f"\n# Calculate the FRIP for {input_fragment_file}")
    snap.metrics.frip(adatas_atac,{"peaks_frac": snap.datasets.cre_HEA()},normalized=True,inplace=True)

    print(f"\n# Calculate the metric summary for each chrom for {input_fragment_file}")
    snap.metrics.summary_by_chrom(adatas_atac, mode='sum')

# --------------------------------------------------------------------------------------------------------------------
#                           FILTERS CELLS BASED ON QC METRICS
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== FILTER CELLS BASED ON QC METRICS =====")
    # Filter cells based on the TSS score and nucleosome signal
    snap.pp.filter_cells(adatas_atac, min_counts=1000, max_counts= 100000,min_tsse=2)
    #adatas_atac = adatas_atac[(adatas_atac.obs['frac_dup'] <= 0.2) & (adatas_atac.obs['frac_mito'] <= 0.3)].copy()

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
    #data.uns['frag_size_distr'] = data.adatas.uns['frag_size_distr'] #calculated at the dataset level
    #data.uns['library_tsse'] = data.adatas.uns['library_tsse']
    #data.uns['frac_overlap_TSS'] = data.adatas.uns['frac_overlap_TSS']
    
    
    

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    data.close()
    
    
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
