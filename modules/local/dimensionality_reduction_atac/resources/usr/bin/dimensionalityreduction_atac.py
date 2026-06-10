#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import argparse                     # command line arguments parser
import re
import warnings
import pathlib                      # library for handle filesystem paths
import numpy as np
import pandas as pd                 # library for data analysis and manipulation
# import scanpy as sc
import snapatac2 as snap
import os
import anndata as ad
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns





warnings.filterwarnings("ignore")
# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function compute the dimensionality reduction for ATAC data using TF-IDF normalization and PCA
    """

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

#Define command line arguments with argparse
    parser = argparse.ArgumentParser(prog='DimensionalityRed', usage='%(prog)s [options]', description = "Dimensionality reduction for ATAC data",
                                    epilog = "This function compute the dimensionality reduction using TF-IDF normalization.",
                                    )
    parser.add_argument('-ad','--input-h5ad-combined',metavar= 'H5AD_INPUT_FILES', type=pathlib.Path, dest='input_h5ad_files',
                        required=True, help="paths of existing matrix files in h5ad format (including file names)")
    parser.add_argument('-fd','--fraction-duplicates', type=float, dest='frac_dup',
                        help="fraction of reads associated with the cell barcode that are duplicates (default: 0.4)", default=0.4)
    parser.add_argument('-pf','--peaks-fraction', type=float, dest='peaks_frac',
                        help="fraction of reads in peaks (default: 0.2)", default=0.2)
    parser.add_argument('-nc','--num-components', type=int, dest='n_comps_atac',
                        help="number of dimensions to keep in dimensionality reduction (default: 30)", default=30)
    parser.add_argument('-nn','--num-neighbours', type=int, dest='n_neighbors_atac',
                        help="number of mutual nearest neighbors used for MNC correction (default: 5)", default=5)
    parser.add_argument('-ncl','--num-clusters', type=int, dest='n_clusters_atac',
                        help="number of clusters used for MNC correction (default: 40)", default=40)
    parser.add_argument('-b', '--blacklist', metavar='BLACKLIST_FILE', type=pathlib.Path, default=None,
                        help="path to the blacklist file in bed format (default is None, no blacklist will be applied)")
    parser.add_argument('-f', '--atac-feature', metavar='N_FEATURES_ATAC', dest='n_features_atac', type=int, default=500000,
                        help="number of most variable features to select for ATAC data (default: 500000)")
    parser.add_argument('-cc', '--cell-counts', metavar='CELL_COUNTS', dest='cell_counts', type=pathlib.Path, default="./cell_counts.csv",
                        help="name of the cell counts csv (default: ./cell_counts.csv)")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.dimred_atac.h5ad",
                        help="path and name of the output h5ad file")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5ad_file = args.input_h5ad_files
    frac_dup = args.frac_dup
    peaks_frac = args.peaks_frac
    n_comps_atac = args.n_comps_atac
    n_neighbors_atac = args.n_neighbors_atac
    n_clusters_atac = args.n_clusters_atac
    blacklist_path = args.blacklist
    n_features_atac = args.n_features_atac
    cell_counts_path = args.cell_counts
    results_dir = args.results
    output =args.out

    # print info on the available matrices
    print("Reading fragment file from the following file:")
    print(f"-File {input_h5ad_file}")


# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print(f"\nProcessing AnnData object in folder {input_h5ad_file} ... ", end ='')

    adata_atac= snap.read(input_h5ad_file,backed=None)
    print("Done!")
    print(f"MuData matrix for combined samples has {adata_atac.shape[0]} cells and {adata_atac.shape[1]} fragments")

# --------------------------------------------------------------------------------------------------------------------
#                           FILTER ANNDATASET OBJECT
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== APPLY FILTERS ON ANNDATASET =====")
    print(f"Cells before any filter: {adata_atac.n_obs}")

    mask = (
        (adata_atac.obs["frac_dup"] < frac_dup) &
        (adata_atac.obs["peaks_frac"] >= peaks_frac)
    )

    print("Cells failing QC filters:", (~mask).sum())

    adata_atac = adata_atac[mask].copy()

    print(f"Cells after all filters: {adata_atac.n_obs}")
    print(adata_atac.obs["sample"].value_counts())

# --------------------------------------------------------------------------------------------------------------------
#                           UPDATE CELL COUNTS
# -------------------------------------------------------------------------------------------------------------------

    # Read cell counts csv to update it
    cell_counts_df = pd.read_csv(cell_counts_path)
    # Save number of cells after second filtering
    counts_after_filters = adata_atac.obs['sample'].value_counts().sort_index()
    cell_counts_df['counts_after_filters'] = cell_counts_df['sample'].map(counts_after_filters).fillna(0).astype(int)

    # Save updated cell counts csv
    cell_counts_df.to_csv(cell_counts_path, index=False)
    print(f"Updated cell counts CSV saved at {cell_counts_path}")

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE FILTERING RESULT
# --------------------------------------------------------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(30,25))
    print("\nVisualize the number of cells after filtering for each sample")
    sns.histplot(adata_atac.obs, x="sample", stat="count", ax=ax)
    locs, labels = plt.xticks()
    ax.set_xlabel("Sample name", fontsize=30)
    ax.set_ylabel("Cell number", fontsize=30)
    plt.setp(labels, rotation=90.,fontsize=30)
    ax.tick_params(axis='y', labelsize=30)
    plt.savefig(os.path.join(results_dir,'Cells_after_filtering_atac.pdf'), bbox_inches='tight', dpi=300)
    plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           DIMENSIONALITY REDUCTION
# --------------------------------------------------------------------------------------------------------------------

    # Compute the spectral embedding for dimensionality reduction
    print("\n===== COMPUTE SPECTRAL EMBEDDING =====")
    snap.pp.add_tile_matrix(adata_atac,bin_size=500,exclude_chroms=["chrM"],min_frag_size=50,max_frag_size=1000,counting_strategy='paired-insertion',inplace=True)


    # Normalize the data using TF-IDF normalization
    print("\nNormalizing data using TF-IDF normalization ... ", end='')
    snap.pp.select_features(adata_atac, n_features=n_features_atac, blacklist=blacklist_path, inplace=True)
    print(adata_atac)
    snap.tl.spectral(adata_atac,n_comps=n_comps_atac,features="selected",weighted_by_sd=True, random_state=0,inplace=True)
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT (before Harmony integration)
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== COMPUTE UMAP BEFORE HARMONY =====")
    print("Computing UMAP before Harmony ... ", end="")
    snap.pp.knn(adata_atac, use_rep="X_spectral")
    snap.tl.umap(adata_atac)
    print("Done!")
        
    print("\n===== PLOT UMAP BEFORE HARMONY =====")
    snap.pl.umap(
        adata_atac,
        color="sample",
        interactive=False,
        show=False,
        out_file=os.path.join(results_dir, "umap_ATAC_sample_before_Harmony.pdf")
    )
    print("Done!")

    # UMAP plot highlighting metadata features (if present)
    # Metadata features have been renamed as meta_* in the convert_mudata step
    pattern = re.compile(r"meta_.*")
    meta_cols = [col for col in adata_atac.obs.columns if pattern.match(col)]
    if len(meta_cols) > 0:
        with PdfPages(os.path.join(args.results, "umap_plot_ATAC_metadata.pdf")) as pdf:
            for col in meta_cols:
                plt.figure(figsize=(45, 35))
                adata_atac.pl.umap(adata_atac, color=col, show=False)
                pdf.savefig(bbox_inches="tight", dpi=300)
                plt.close()
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           BATCH CORRECTION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== PERFORM BATCH CORRECTION =====")
    print("Performing batch correction with Harmony... ", end='')
    X_harmony = snap.pp.harmony(
        adata_atac,
        batch="sample",
        use_rep="X_spectral",
        inplace=False,
    )
    X_harmony = np.asarray(X_harmony)
    if X_harmony.shape[0] != adata_atac.n_obs and X_harmony.shape[1] == adata_atac.n_obs:
        X_harmony = X_harmony.T
    assert X_harmony.shape[0] == adata_atac.n_obs, X_harmony.shape
    adata_atac.obsm["X_spectral_harmony"] = X_harmony
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           CLUSTERING ON TILE MATRIX
# --------------------------------------------------------------------------------------------------------------------

    # Perform clustering
    print("\n===== PERFORM CLUSTERING =====")
    print(f"Performing clustering ... \n", end='')
    snap.pp.knn(adata_atac, use_rep="X_spectral_harmony")
    resolutions = np.round(np.arange(0.1, 1.1, 0.1), 2)
    clustering_labels = []
    for res in resolutions:
        clustering_labels.append("leiden_{}".format(res))
        print("calculating leiden_{}".format(res))
        snap.tl.leiden(adata_atac, resolution=res, key_added="leiden_tile_{}".format(res))
    print("Done!")

    print(adata_atac)
    
# --------------------------------------------------------------------------------------------------------------------
#                           COMPUTE AND VISUALIZE UMAP PLOT (after Harmony integration)
# --------------------------------------------------------------------------------------------------------------------

    # Compute UMAP for visualization
    print("\n===== COMPUTE UMAP =====")
    print("Computing UMAP ... ", end='')
    snap.tl.umap(adata_atac, min_dist=0, use_rep="X_spectral_harmony")

    print("Done!")

    print("\nVisualize UMAP plots after integration")
    snap.pl.umap( adata_atac, color="sample", interactive=False, show=False, out_file=os.path.join(results_dir, "umap_ATAC_sample_Harmony.pdf"))
    
    # Visualize UMAP for all Leiden resolutions
    print("\n===== VISUALIZE UMAP FOR ALL RESOLUTIONS =====")
    for res in resolutions:
        leiden_key = f"leiden_tile_{res}"
        print(f"Visualizing {leiden_key}")
        print(f"\n      Visualizing UMAP colored by {leiden_key}")
        snap.pl.umap(
            adata_atac,
            color=leiden_key,
            interactive=False,
            show=False,
            out_file=os.path.join(results_dir, f"umap_ATAC_res_{res}_Harmony.png")
        )
    print("Done!")

    print("\n===== PRINT ALL UMAPS TOGETHER =====")
    image_files = [
        os.path.join(results_dir, f"umap_ATAC_res_{res}_Harmony.png")
        for res in resolutions
    ]
    n_images = len(image_files)
    n_cols = 4
    n_rows = int(np.ceil(n_images / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(5 * n_cols, 5 * n_rows),
        squeeze=False
    )
    axes = axes.flatten()
    for i, img_path in enumerate(image_files):
        ax = axes[i]
        img = plt.imread(img_path)
        ax.imshow(img)
        ax.axis("off")
    # Remove unused axes
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout()
    plt.savefig(
        os.path.join(results_dir, "umap_ATAC_all_res_Harmony.pdf"),
        bbox_inches="tight",
        dpi=300
    )
    plt.close()
    print("Done!")
    
    print("\n===== PRINT UMAPS BY METADATA =====")
    # UMAP plot highlighting metadata features (if present)
    # Metadata features have been renamed as meta_* in the convert_mudata step
    pattern = re.compile(r"^meta_.*")
    meta_cols = [col for col in adata_atac.obs.columns if pattern.match(col)]
    print(f"Saving umaps for metadata{meta_cols}")
    if len(meta_cols) > 0:
        with PdfPages(os.path.join(results_dir, "umap_ATAC_meta_Harmony.pdf")) as pdf:
            for col in meta_cols:
                if col in adata_atac.obs.columns:
                    print(f"Visualizing UMAP colored by {col}")
                    snap.pl.umap(
                        adata_atac,
                        color=col,
                        interactive=False,
                        show=False,
                        out_file=os.path.join(results_dir, f"umap_ATAC_meta_{col}_Harmony.pdf")
                    )
                    plt.close()
    print("Done!")
    
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING OUTPUT FILE =====")
    print(f"Saving atac data in {output}")
    adata_atac.write(output)
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
