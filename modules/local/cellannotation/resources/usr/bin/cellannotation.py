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
    parser.add_argument('-csv_annotation', '--csv_annotation_out', metavar='METADATA', default="metadata.csv",
                        help="path and name of csv table with metadata information for each cell")
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
    output_csv_annotation= args.csv_annotation_out
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
    gex.var['gene_name'] = gex.var['gene_name'].astype(str)
    gex.var['gene_id'] = gex.var.index.astype(str)
    gex.var = gex.var.set_index('gene_name')
    print(gex.var)



# -------------------------------------------------------------------------------------------------------------------
#                                 CELLTYPIST ANNOTATION
# --------------------------------------------------------------------------------------------------------------------
    with open(input_model_list, 'r') as f:
        model_names = [line.strip() for line in f if line.strip()]

    df_list = []

    for model_name in model_names:
        print(f"\n=== Annotating with model: {model_name} ===")
        model = load_model(model_name)


        predictions = celltypist.annotate(gex, model=model, majority_voting=True,mode = 'best match')
        predictions_adata = predictions.to_adata()

        model_base = model_name.replace(".pkl", "").replace(" ", "_")
        col_name = f"{model_base}:majority_voting"
        col_conf = f"{model_base}:conf"

        df_celltypist = predictions_adata.obs.loc[
            gex.obs.index, ["majority_voting", "conf_score"]
        ].rename(columns={
            "majority_voting": col_name,
            "conf_score": col_conf
        })

        df_celltypist.columns = [f"celltypist:{col_name}", f"celltypist:{col_conf}"]
        df_list.append(df_celltypist)

    df_celltypist = pd.concat(df_list, axis=1)
    print(df_celltypist)
    gex.obs = pd.concat([gex.obs, df_celltypist], axis=1)
    cols_to_drop = [
            'predicted_labels',   # raw cell-level predicted label before majority voting
            'majority_voting',    # label after local subcluster majority voting (already renamed and saved)
            'conf_score',         # confidence score (already renamed and saved)
            'over_clustering'     # internal over-clustering information used by CellTypist
    ]

    gex.obs = gex.obs.drop(columns=[c for c in cols_to_drop if c in gex.obs.columns])
    print(gex.obs)

# --------------------------------------------------------------------------------------------------------------------
#                           SUMMARY OF CELLTYPIST ANNOTATION
# --------------------------------------------------------------------------------------------------------------------

    for model_name in model_names:
        clean_model_name = model_name.replace(".pkl", "").replace(" ", "_")
        col_name = f"celltypist:{clean_model_name}:majority_voting"

        summary_table_pool = gex.obs.groupby(['sample', col_name]).size().reset_index(name='count')

        output_csv_pool = output_csv.with_name(f"{output_csv.stem}_{clean_model_name}_by_pool.csv")
        summary_table_pool.to_csv(output_csv_pool, index=False)

        # Per Inferred_donor (idem)
        #summary_table_sample = gex.obs.groupby(['Inferred_donor', col_name]).size().reset_index(name='count')
        #output_csv_sample = output_csv.with_name(f"{output_csv.stem}_{clean_model_name}_by_sample.csv")
        #summary_table_sample.to_csv(output_csv_sample, index=False)

        print(f"Saved summaries for model {clean_model_name}")

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------
        print("\n===== VISUALIZING UMAP PLOT =====")
        print(f"\nVisualized batch-corrected UMAP plot for model {clean_model_name}")
        # Visualize batch-corrected UMAP plot
        plt.figure(figsize=(35, 25))
        mu.pl.embedding(gex, color =col_name, basis= 'X_umap', show=False)
        plt.savefig(os.path.join(args.results, f"Annotated_UMAP_{clean_model_name}.pdf"), bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Saved summaries and UMAP for model {clean_model_name}")


# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
    gex.var.reset_index(inplace=True)
    gex.var.index.name = None
    mdata.mod['gex'] = gex
    mdata.update()
    print(mdata.obs)

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE METADATA INFORMATION INTO A EXCEL FILE
# --------------------------------------------------------------------------------------------------------------------

    mdata.obs.to_csv(output_csv_annotation, index=True)
    print("Metadata information for each cell saved in csv file {}".format(output_csv_annotation))


# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print(mdata.obs)
    print(mdata.var)

    print("Saving h5ad data to file {}".format(output))
    mdata.write(output)
    print(mdata)


#####################################################################################################


if __name__ == '__main__':
    main()
