#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import os                           # filesystem utilities
import pathlib                      # library for handle filesystem paths
import matplotlib.pyplot as plt     # library for visualization
import pandas as pd                 # library for data analysis and manipulation
import scanpy as sc                 # single-cell data processing
import decoupler as dc
import pertpy as pt
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
    This function performs pseudobulks and differential expression (DE) analysis.
    """
# --------------------------------------------------------------------------------------------------------------------
#                                          LIBRARY CONFIG
# --------------------------------------------------------------------------------------------------------------------

    sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

# Define command line arguments with argparse

    parser = argparse.ArgumentParser(prog='differential_exp', usage='%(prog)s [options]',description = "Pseudobulk and differential expression analysis",
        epilog = "This function computes pseudobulk profiles within each group and sample of interest and computes the DGE analysis.")
    parser.add_argument('-ad','--input-h5mu-file',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    #parser.add_argument('-c','--csv_comparison',metavar= 'COMPARISON_INPUT_CSV', type=pathlib.Path, dest='comparison_csv',
    #                    required=True, help="paths of file of comparison groups in csv format")
    #parser.add_argument('-f','--csv_formula',metavar= 'COMPARISON_INPUT_FORMULA', type=pathlib.Path, dest='formula_csv',
    #                    required=True, help="paths of file of formula in csv format")
    parser.add_argument('-e', '--excel_out', metavar='DIFFERENTIAL_GENES_XLSX', default="differential_genes.xlsx",
                        help="path and name of excel table with differential genes")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_file = args.input_h5mu_files
    output_xlsx=args.excel_out
    #comparison_csv=args.comparison_csv
    #formula_csv=args.formula_csv

    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print(f"-File {input_h5mu_file}:")

    #print("Reading comparison groups from the following file:")
    #print(f"-File {comparison_csv}:")

    #print("Reading formula from the following file:")
    #print(f"-File {formula_csv}:")
# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the MTX combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder ... ", end ='')
    mdata= md.read(input_h5mu_file)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ COMPARISON AND FORMULA FILE
# --------------------------------------------------------------------------------------------------------------------

    #comparison = pd.read_csv(comparison_csv).to_string(index=False)
    #print(comparison)


    #with open(formula_csv, 'r') as file:
    # Read the first line and store it as a string
    #    content = str(file.readline())
# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    
    print("\n===== GEX MODALITY DATA =====")
    gex = mdata.mod['gex']
    gex.var['gene_symbols'] = gex.var['gene_symbols'].astype('str')
    gex.var_names = gex.var["gene_symbols"]
    gex.var_names_make_unique()
    
    
    amplified_samples = ["SK-N-DZ", "CHP134","BE2c"]
    no_amplified_samples = ["SHSY5Y", "SK-N-SH", "SK-N-AS"]

    gex.obs["MYCN_status"] = "unknown"
    gex.obs.loc[gex.obs['sample'].isin(amplified_samples), "MYCN_status"] = "amplified"
    gex.obs.loc[gex.obs['sample'].isin(no_amplified_samples), "MYCN_status"] = "no_amplified"



# --------------------------------------------------------------------------------------------------------------------
#                                 PSEUDOBULK
# --------------------------------------------------------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(40,10))
    print("\n===== PSEUDOBULK CREATION =====")
    # extract pseudobulk profiles within each group and sample of interest
    print("\nExtracting the pseudobulk profile for each sample and condition")
    ps = pt.tl.PseudobulkSpace()
    pdata = ps.compute(gex, target_col="sample", groups_col="MYCN_status", layer_key="count", mode="sum", min_cells=10, min_counts=1000)
    print(pdata)
    ps.plot_psbulk_samples(pdata, groupby=["sample", "MYCN_status"], figsize=(12, 4),return_fig=True)
    plt.savefig(os.path.join(args.results,'Pseudobulk_counts.png'))
    plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           DIFFERENTIAL EXPRESSION ANALYSIS
# --------------------------------------------------------------------------------------------------------------------
 

    #formula_str = '~Tissue'
    #print(type(content))
    pds2 = pt.tl.PyDESeq2(pdata, design='MYCN_status')
    print("\nFitting the model and testing the contrasts")
    pds2.fit()
    print(pds2.design.head())

  
    res_df = pds2.test_contrasts(pds2.contrast(column='MYCN_status', baseline="no_amplified", group_to_compare="amplified"))
    res_df.head(10)
    
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE DIFFERENTIAL GENES INTO EXCEL FILE
# --------------------------------------------------------------------------------------------------------------------

    print("\nSaving differential genes into excel file")
    res_df.to_excel(output_xlsx, index=False)
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE COMPARISON RESULTS
# --------------------------------------------------------------------------------------------------------------------

    res_df = pds2.compare_groups(pdata, column="MYCN_status", baseline="no_amplified", groups_to_compare="amplified")
    #edgr.plot_multicomparison_fc(res_df, figsize=(12, 1.5))
    print(res_df)

    # Visualize volcano plot
    fig, ax = plt.subplots(figsize=(40,10))
    print("\nVisualized volcano plot")
    pds2.plot_volcano(res_df, log2fc_thresh=0,return_fig=True,figsize=(10, 4))
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90.)
    plt.savefig(os.path.join(args.results,'Volcano_plot.png'))
    plt.close()

    print("\nVisualized fold changes of the top differentially expressed genes")
    pds2.plot_fold_change(res_df, n_top_vars=15,return_fig=True,figsize=(14, 8))
    plt.savefig(os.path.join(args.results,'Fold_change.png'))
    plt.close()

#####################################################################################################
    

if __name__ == '__main__':
    main()
