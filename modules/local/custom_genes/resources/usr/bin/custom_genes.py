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
import seaborn as sns               # library for statistical data visualization
import pandas as pd                 # library for data analysis and manipulation
import scanpy as sc                 # single-cell data processing
import scanpy.external as sce       # library for harmony integration
import muon as mu
import numpy as np
import anndata as ad
from matplotlib.backends.backend_pdf import PdfPages

warnings.filterwarnings("ignore")
# PARAMETERS

# set script version number
VERSION = "0.0.1"

# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function clusters single-cell data at different resolution levels.
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
    parser = argparse.ArgumentParser(description = "Custom genes plotting",
        epilog = "This function produces plots to highlight the expression of a list of custom genes exernally provuded.",
        )
    parser.add_argument('-ad','--input-h5mu-file',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('-g', '--input-genelist', metavar='INPUT_CUSTOM_GENES', default="",
                        help="path of a .txt file containing a list of custom genes to plot (one gene per line, no header)")
    parser.add_argument('-res', '--resolution',  dest='set_res', type=float, default=100,
                        help="clustering resolution. By default, no resolution is considered.")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()


# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_file = args.input_h5mu_files
    input_genelist = args.input_genelist
    set_res = args.set_res

    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print("-File {}:".format(str(input_h5mu_file)))

# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5MU FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the MTX combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder ... ", end ='')
    mdata = mu.read_h5mu(input_h5mu_file)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    # Extracting gex modality data
    print("\n===== GEX MODALITY DATA =====")
    gex = mdata.mod['gex']

# --------------------------------------------------------------------------------------------------------------------
#                                 INPUT GENE LIST
# --------------------------------------------------------------------------------------------------------------------
    # Importing external gene list
    print("\n===== READING GENE LIST =====")
    geneset = pd.read_csv(input_genelist, header=None)[0].tolist()
    print(f"Imported genelist with {len(geneset)} genes.")
    geneset = [g for g in geneset if g in gex.var_names]
    print(f"Among these, {len(geneset)} genes are expressed and will be used for plotting.")
    # Name to be used for saving output files
    geneset_name = os.path.splitext(os.path.basename(input_genelist))[0]
    # Check if any gene from the list is expressed
    if len(geneset) == 0:
        print("No genes from the provided list are expressed in the dataset. Exiting.")
        return

# --------------------------------------------------------------------------------------------------------------------
#                                 FEATURES PLOTS
# --------------------------------------------------------------------------------------------------------------------
    # Scoring gene set
    sc.tl.score_genes(gex, gene_list=geneset, score_name="geneset_score")

    # Feature plotting for custom gene set, both single genes and overall geneset score
    print("\n===== FEATURE PLOTS =====")
    feat_pdf = f"{geneset_name}_features_plots.pdf"
    with PdfPages(os.path.join(args.results, feat_pdf)) as pdf_feature:
        for gene in geneset:
            plt.figure(figsize=(45, 35))
            sc.pl.umap(gex, color=gene, title=f"Expression of {gene}", show=False)
            pdf_feature.savefig()
            plt.close()
        plt.figure(figsize=(45, 35))
        sc.pl.umap(gex, color="geneset_score", title=f"{geneset_name} overall score", show=False)
        pdf_feature.savefig()
        plt.close()

    print("\Features plots completed successfully.")

# --------------------------------------------------------------------------------------------------------------------
#                              VIOLIN, DOTPLOT AND HEATMAP WITH SPECIFIED RESOLUTION
# --------------------------------------------------------------------------------------------------------------------
    # Producing dotplot and heatmap for custom gene set at specified resolution
    print("\n===== HEATMAP, VIOLIN AND DOTPLOT =====")
    if set_res != 100:
        cluster_key = "leiden_{}".format(set_res)
        if cluster_key not in gex.obs.keys():
            print("Clustering at resolution {} not found. Exiting.".format(set_res))
            return

        # Dotplot
        dotplot_pdf = f"{geneset_name}_dotplot_r{set_res}.pdf"
        with PdfPages(os.path.join(args.results, dotplot_pdf)) as pdf_dotplot:
            plt.figure(figsize=(55, 35))
            sc.pl.dotplot(gex, var_names=geneset, groupby=cluster_key, show=False, standard_scale='var', title=f"{geneset_name}")
            pdf_dotplot.savefig()
            plt.close()
        print("\Dotplot completed successfully.")

        # Heatmap
        heatmap_pdf = f"{geneset_name}_heatmap_r{set_res}.pdf"
        with PdfPages(os.path.join(args.results, heatmap_pdf)) as pdf_heatmap:
            plt.figure(figsize=(45, 55))
            sc.pl.heatmap(gex, var_names=geneset, groupby=cluster_key, show_gene_labels=True, show=False)
            plt.title(f"{geneset_name}")
            pdf_heatmap.savefig()
            plt.close()
        print("\Heatmap completed successfully.")

        # Violin plot
        violin_pdf = f"{geneset_name}_violin_r{set_res}.pdf"
        with PdfPages(os.path.join(args.results, violin_pdf)) as pdf_violin:
            for gene in geneset:
                plt.figure(figsize=(45, 35))
                sc.pl.violin(gex, keys=gene, groupby=cluster_key, show=False)
                plt.title(f"Expression of {gene}")
                pdf_violin.savefig()
                plt.close()
            plt.figure(figsize=(45, 35))
            sc.pl.violin(gex, keys="geneset_score", groupby=cluster_key, show=False)
            plt.title(f"{geneset_name} overall score")
            pdf_violin.savefig()
            plt.close()
        print("\Violin plot completed successfully.")



#####################################################################################################

# Actual execution
if __name__ == '__main__':
    main()
