#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================


# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import os                           # filesystem utilities
import re                           # hanlding regex
import pathlib                      # library for handle filesystem paths
import matplotlib.pyplot as plt     # library for visualization
import pandas as pd                 # library for data analysis and manipulation
import scanpy as sc                 # single-cell data processing
import scanpy.external as sce       # library for harmony integration
import mudata as md
import muon as mu
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
    This function integrates single-cell data from multiple experiments, using the selected integration variable (e.g. 'sample','patient', etc.).
    If the integration variable contains only one unique level (e.g. a single sample),
    Harmony integration is skipped and the existing UMAP from the dimensionality
    reduction step is reused. All output files are produced regardless.
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


    parser = argparse.ArgumentParser(prog='Int', usage='%(prog)s [options]', description = "Data integration",
        epilog = "This function integrate the single-cell dataset based on run_id.")
    parser.add_argument('-ad', '--input-h5mu-file', metavar='H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5MU_OUTPUT_FILE', type=pathlib.Path, default="matrix.integrated.h5mu",
                        help="name of the output h5ad file after integration")
    parser.add_argument('-csv', '--csv_out', metavar='CSV_TABLE', type=pathlib.Path, default="Harmony_UMAP_coordinates_GEX.csv",
                        help="path and name of csv table with UMAP coordinates for each cell")
    parser.add_argument('-nnh', '--n_neighbors_harmony', dest='n_neighbors_harmony', type=int, default=20, help="Size of local neighborhood used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. Values should be in the range 2 to 100")
    parser.add_argument('-mdh', '--min_dist_harmony', dest='min_dist_harmony', type=float, default=0.1, help="minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together")
    parser.add_argument('-var','--integration_var', type=str, default='sample', help="variable in .obs to use for integration (default is 'sample')")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()


# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------


    print("\n===== INPUT H5AD FILES =====")
    input_h5mu_file = args.input_h5mu_files
    output = args.out
    output_csv= args.csv_out
    n_neighbors_harmony = args.n_neighbors_harmony
    min_dist_harmony = args.min_dist_harmony
    integration_var = args.integration_var


    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print(f"-File {input_h5mu_file}:")


# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES
# --------------------------------------------------------------------------------------------------------------------


    # Read folders with the MTX combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder ... ", end ='')
    mdata = md.read(input_h5mu_file)
    print("Done!")
    print(f"Count matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")


# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== GEX MODALITY DATA =====")
    gex = mdata.mod['gex']


# --------------------------------------------------------------------------------------------------------------------
#                                 DATA INTEGRATION
# --------------------------------------------------------------------------------------------------------------------


    print("\n===== DATA INTEGRATION =====")
    # Check if integration_var is in obs; try meta_ prefix as fallback
    if integration_var not in gex.obs.columns:
        if f"meta_{integration_var}" in gex.obs.columns:
            print(f"Variable '{integration_var}' not found in .obs. Using 'meta_{integration_var}' instead.")
            integration_var = f"meta_{integration_var}"
        else:
            print(f"Error: Variable '{integration_var}' not found in .obs columns.")
            print(f"Available columns: {list(gex.obs.columns)}")
            raise ValueError(f"Integration variable '{integration_var}' not found in .obs")

    # Check how many unique levels the integration variable has
    n_levels = gex.obs[integration_var].nunique()
    print(f"Integration variable '{integration_var}' has {n_levels} unique level(s): "
          f"{sorted(gex.obs[integration_var].unique().tolist())}")

    if n_levels > 1:
        # ------------------------------------------------------------------
        # MULTI-LEVEL: run Harmony integration and recompute neighbors/UMAP
        # ------------------------------------------------------------------
        print("\nMultiple levels detected — performing Harmony integration.")
        sce.pp.harmony_integrate(gex, integration_var)

        print("\n===== BATCH-CORRECTED UMAP =====")
        print("\nConstruction of the nearest neighbor graph (Harmony embedding)")
        sc.pp.neighbors(gex, n_neighbors=n_neighbors_harmony, use_rep="X_pca_harmony")

        print("\nPerforming UMAP on Harmony-corrected embedding")
        sc.tl.umap(gex, min_dist=min_dist_harmony, random_state=42)

    else:
        # ------------------------------------------------------------------
        # SINGLE-LEVEL: skip Harmony, reuse the UMAP from the dim-red step
        # ------------------------------------------------------------------
        print(f"\nOnly one unique level found in '{integration_var}' — "
              "skipping Harmony integration.")
        print("Reusing UMAP coordinates computed during dimensionality reduction.")

        # Sanity check: make sure the UMAP was already computed upstream
        if "X_umap" not in gex.obsm:
            raise RuntimeError(
                "X_umap not found in gex.obsm. "
                "Run the dimensionality reduction module before the integration module."
            )


# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------


    # Visualize batch-corrected UMAP plot (or passthrough UMAP for single sample)
    print("\nVisualized batch-corrected UMAP plot")
    plt.figure(figsize=(45,35))
    mu.pl.embedding(gex, color ='sample', basis= 'X_umap', show=False)
    plt.savefig(os.path.join(args.results,'Harmony_corrected_UMAP_plot_GEX.pdf'), bbox_inches='tight', dpi=300)
    plt.close()


    # UMAP plots highlighting n_genes and mito_percent distribution
    plt.figure(figsize=(55, 55))
    sc.pl.umap(gex, color=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"], wspace=0.5, ncols=2)
    plt.savefig(os.path.join(args.results,'Harmony_corrected_UMAP_plot_GEX_QC.pdf'), bbox_inches='tight', dpi=300)
    plt.close()


    # UMAP plot highlighting cell cycle phases (if present)
    if 'phase' in gex.obs.columns:
        plt.figure(figsize=(45, 35))
        mu.pl.embedding(gex, color='phase', basis='X_umap', show=False)
        plt.savefig(os.path.join(args.results,'Harmony_corrected_UMAP_plot_GEX_phase.pdf'), bbox_inches='tight', dpi=300)
        plt.close()


    # UMAP plot highlighting celltypist annotation (if present)
    pattern = re.compile(r"celltypist:.*:majority_voting")
    celltypist_cols = [col for col in gex.obs.columns if pattern.match(col)]
    if len(celltypist_cols) > 0:
        with PdfPages(os.path.join(args.results, "Harmony_corrected_UMAP_plot_GEX_celltypist.pdf")) as pdf:
            for col in celltypist_cols:
                plt.figure(figsize=(45, 35))
                mu.pl.embedding(gex, color=col, basis="X_umap", show=False)
                pdf.savefig(bbox_inches="tight", dpi=300)
                plt.close()


    # UMAP plot highlighting metadata features (if present)
    # Metadata features have been renamed as meta_* in the convert_mudata step
    pattern = re.compile(r"meta_.*")
    meta_cols = [col for col in gex.obs.columns if pattern.match(col)]
    if len(meta_cols) > 0:
        with PdfPages(os.path.join(args.results, "Harmony_corrected_UMAP_plot_GEX_metadata.pdf")) as pdf:
            for col in meta_cols:
                plt.figure(figsize=(45, 35))
                mu.pl.embedding(gex, color=col, basis="X_umap", show=False)
                pdf.savefig(bbox_inches="tight", dpi=300)
                plt.close()


# --------------------------------------------------------------------------------------------------------------------
#                           SAVE GEX DATA INTO MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")
    mdata.mod['gex'] = gex
    mdata.update()


# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    print(f"Saving h5mu data to file {output}")
    mdata.write(output)
    print("Done!")


    df = pd.DataFrame(gex.obsm["X_umap"], index=gex.obs_names).rename(columns={0: "X_UMAP", 1: "Y_UMAP"})
    df.index.name = 'cell_barcodes'
    print(f"Saving csv table with Harmony corrected UMAP coordinates for each cell {output_csv}")
    df.to_csv(output_csv)
    print("Done!")


#####################################################################################################


if __name__ == '__main__':
    main()
