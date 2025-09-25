#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT
import os
import platform
import warnings
import argparse
import pathlib
from pathlib import Path
import pandas as pd
import scanpy as sc
import muon as  mu
import matplotlib.pyplot as plt
import celltypist
from celltypist import models as ct_models
import pickle

os.environ["CELLTYPIST_OFFLINE"] = "1"

warnings.filterwarnings("ignore")
# PARAMETERS

# set script version number
VERSION = "0.0.1"

# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================
def load_model(model_name):
    model_file = f"{model_name}.pkl" if not str(model_name).endswith(".pkl") else model_name

    # Check if the model file exists in the current directory
    if os.path.exists(model_file):
        print(f"Loading model {model_file} from local directory.")
        model_obj = ct_models.Model.load(str(model_file))
    else:
        # If the model file does not exist, download it
        print(f"Model {model_file} not found locally, downloading...")
        ct_models.download_models(model=model_file)
        model_obj = ct_models.Model.load(str(model_file))

    return model_obj


def main():
    """
    This function annotates scRNA-seq query data by retrieving the most likely cell type labels from either the built-in CellTypist models.
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

    parser = argparse.ArgumentParser(description = "CellTypist",
        epilog = "This function annotates scRNA-seq query data by retrieving the most likely cell type labels from either the built-in CellTypist models",
        )
    parser.add_argument('-ad','--input-h5mu-file',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('--model-list', metavar='MODEL_LIST_TXT', type=pathlib.Path, required=True,
                    help="path to a .txt file listing the CellTypist model .pkl files to use (one per line)")
    parser.add_argument('-o', '--out', metavar='H5MU_OUTPUT_FILE', type=pathlib.Path, default="matrix.annotated.h5mu",
                        help="name of the output h5ad file after cell annotation")
    parser.add_argument('-csv', '--csv_out', metavar='CELL_ANNOTATION', default="summary_cellannotation.csv",
                        help="path and name of csv table cell annotation summary")
    parser.add_argument('-e', '--excel_out', metavar='METADATA', default="metadata.xlsx",
                        help="path and name of excel table with metadata information for each cell")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_file = args.input_h5mu_files
    input_model_list = args.model_list
    output_csv= Path(args.csv_out)
    output_excel= args.excel_out
    output = args.out


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
    mdata= mu.read_h5mu(input_h5mu_file)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== GEX MODALITY DATA =====")
    gex = mdata.mod['gex']
    gex.var['gene_symbols'] = gex.var['gene_symbols'].astype(str)
    gex.var = gex.var.set_index('gene_symbols')
    print(gex.var)



# -------------------------------------------------------------------------------------------------------------------
#                                 CELLTYPIST ANNOTATION
# --------------------------------------------------------------------------------------------------------------------

    df_list = []

    model = load_model(input_model_list)
    model_name = input_model_list.stem

    predictions = celltypist.annotate(gex, model=model, majority_voting=True,mode = 'prob match', p_thres = 0.5)
    predictions_adata = predictions.to_adata()

    df_celltypist = predictions_adata.obs.loc[
        gex.obs.index, ["predicted_labels", "conf_score"]
    ]

    df_celltypist.columns = [f"celltypist:{model_name}", f"celltypist:{model_name}:conf"]
    df_list.append(df_celltypist)

    df_celltypist = pd.concat(df_list, axis=1)

    gex.obs = pd.concat([gex.obs, df_celltypist], axis=1)

# --------------------------------------------------------------------------------------------------------------------
#                           SUMMARY OF CELLTYPIST ANNOTATION
# --------------------------------------------------------------------------------------------------------------------
    output_csv_pool = output_csv.with_name(output_csv.stem + "_by_pool.csv")
    summary_table_pool = gex.obs.groupby(['sample','predicted_labels']).size().reset_index(name='count')
    print(summary_table_pool)
    summary_table_pool.to_csv(output_csv_pool,index=False)
    print("Done!")

    output_csv_sample = output_csv.with_name(output_csv.stem + "_by_sample.csv")
    summary_table_sample = gex.obs.groupby(['Inferred_donor', 'predicted_labels']).size().reset_index(name='count')
    print(summary_table_sample)
    summary_table_sample.to_csv(output_csv_sample,index=False)
    print("Done!")
# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

    # Visualize batch-corrected UMAP plot
    print("\nVisualized batch-corrected UMAP plot")
    mu.pl.embedding(gex, color =['predicted_labels'],basis= 'X_umap',legend_loc='on data',show=False)
    plt.savefig(os.path.join(args.results,'Annotated_UMAP_plot_GEX.png'))
    plt.close()


# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
    mdata.mod['gex'] = gex
    mdata.update()
    print(mdata.obs)

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE METADATA INFORMATION INTO A EXCEL FILE
# --------------------------------------------------------------------------------------------------------------------

    mdata.obs.to_excel(output_excel, index=False)
    print("Metadata information for each cell saved in excel file {}".format(output_excel))


# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print(mdata.obs)
    print("Saving h5ad data to file {}".format(output))
    mdata.write(output)
    print(mdata)


#####################################################################################################


if __name__ == '__main__':
    main()
