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
        sorted_by_barcode=False
    )
    print(adatas_atac)
    
    #adata_atac = snap.pp.import_fragments(input_fragment_file,file=output,chrom_sizes=snap.genome.GRCh38,sorted_by_barcode=False)
    
    #print(f"MuData matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")
    
# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE AND VISUALIZE QUALITY METRICS
# --------------------------------------------------------------------------------------------------------------------
        
    # Compute the fragment size distribution for each sample

    print("\n===== COMPUTE QUALITY METRICS {} =====")
    print(f"\n# Calculate the fragment size distribution for {input_fragment_file}")

    #fig = snap.pl.frag_size_distr(adatas_atac, show=False)
    #fig.update_yaxes(type="log")
    #fig.write_image("Nucleosome_signal_1.png", width=1000, height=600)

    print(f"\n# Calculate the TSS score for {input_fragment_file}")
    # Compute the TSS score for each sample
    #snap.metrics.tsse(adatas_atac, snap.genome.hg38)
    #fig2 = snap.pl.tsse(adata_atac, interactive=False)
    #fig2.write_image("TSS_score_1.png", width=1000, height=600)

    print(f"\n# Calculate the FRIP for {input_fragment_file}")
    #snap.metrics.frip(adatas_atac,{"peaks_frac": snap.datasets.cre_HEA()},inplace=True)

    snap.pp.add_tile_matrix(adatas_atac,bin_size=500)
    snap.pp.select_features(adatas_atac, n_features=100000)
    snap.pp.scrublet(adatas_atac)
    
    
    
    data = snap.AnnDataSet(
    adatas=[(name, adata) for (name, _), adata in zip(files, adatas_atac)],
    filename=output
    )
    print(data)
    
    
    unique_cell_ids = [sa + ':' + bc for sa, bc in zip(data.obs['sample'], data.obs_names)]
    data.obs_names = unique_cell_ids
    assert data.n_obs == np.unique(data.obs_names).size
    

    #adata_merged = ad.concat(adatas_atac, label="sample", keys=[name for name, _ in files])
    #adata_merged.obs_names = [f"{name}_{bc}" for name, bc in zip(adata_merged.obs["sample"], adata_merged.obs_names)]


    
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    data.close()
    #adata_merged.write(output)
    
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
