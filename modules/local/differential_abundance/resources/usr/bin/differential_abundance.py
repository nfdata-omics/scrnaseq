#!/usr/bin/env python

###############################################################################
######################## Differential abundance analysis #######################
###############################################################################
# Following vignette: https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/milo.html
# Import libraries
import warnings
import argparse
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pertpy as pt
import mudata as md
from matplotlib.backends.backend_pdf import PdfPages

warnings.filterwarnings("ignore")

# Set to not interactive
matplotlib.use("Agg")

#########################################
# Define and parse CLI arguments
#########################################
parser = argparse.ArgumentParser(description="Differential abundance analysis")

parser.add_argument("--mdata", type=str, required=True,
                    help="Path to input MuData object (h5mu file)")
parser.add_argument("--target", type=str, required=True,
                    help="Name of the target condition (e.g. treated)")
parser.add_argument("--reference", type=str, required=True,
                    help="Name of the reference condition (e.g. control)")
parser.add_argument("--column_to_test", type=str, required=True,
                    help="Name of the column to use for differential abundance (e.g. condition)")

args = parser.parse_args()


#########################################
# Load data
#########################################
# Set workdir
os.makedirs('diff_abundance', exist_ok=True)

# Import object
print("Loading MuData object...")
mdata = md.read(args.mdata)
gex = mdata.mod['gex']

# Handle gene names (parse version)
gex.var['gene_name'] = gex.var['gene_name'].astype('str')
gex.var_names = gex.var["gene_name"]
gex.var_names_make_unique()

# Set categories to compare
target_level = args.target
ref_level = args.reference

# Set real name of column with target and reference
column_name = "meta_" + args.column_to_test

# Check how many samples per group
gex.obs[["sample", column_name]].drop_duplicates().value_counts(column_name)

# Create Excel file with how many cells per sample per cluster / per group per cluster
with pd.ExcelWriter(os.path.join(os.getcwd(),'diff_abundance/cell_counts.xlsx')) as writer:
    cell_counts = gex.obs.groupby(["sample", "leiden_0.5"]).size().unstack(fill_value=0)
    cell_counts.to_excel(writer, sheet_name="samples_vs_clusters")
    cell_counts = gex.obs.groupby([column_name, "leiden_0.5"]).size().unstack(fill_value=0)
    cell_counts.to_excel(writer, sheet_name="group_vs_clusters")

# Initialize object for Milo analysis: mudata object storing rna e milo matrix for differential abundance
milo = pt.tl.Milo()
mmdata = milo.load(gex)
mmdata['milo'].obs['sample'] = mmdata['rna'].obs['sample']
mmdata['milo'].obs[column_name] = mmdata['rna'].obs[column_name]
mmdata['milo'].obs['leiden_0.5'] = mmdata['rna'].obs['leiden_0.5']

# Check if knn graph is already compiled (it should be, since the object has been already processed with Scanpy workflow)
print("KNN graph already compiled:", "neighbors" in gex.uns) # True, ok

#########################################
# Load data
#########################################
# Construct neighbourhoods
# Group cells into clusters or local regions on the KNN graph, which will be used to analyze differences between samples or conditions
# Use 10% of all cells as centres to build neighbourhoods
milo.make_nhoods(mmdata["rna"], prop=0.1)
# number of neighbourhoods
mmdata["rna"].obsm["nhoods"]
# The information on which cells are sampled as index cells of representative neighbourhoods is stored in mdata['rna'].obs, along with the distance of the index to the kth nearest neighbor, which is used later for the SpatialFDR correction.
mmdata["rna"][mmdata["rna"].obs["nhood_ixs_refined"] != 0].obs[["nhood_ixs_refined", "nhood_kth_distance"]]

# Plot the distribution of neighbourhood sizes, to check that the minimal value of k makes sense, and that the distribution of sizes is not too wide
nhood_size = np.array(mmdata["rna"].obsm["nhoods"].sum(0)).ravel()
plt.figure(figsize=(6, 4))
plt.hist(nhood_size, bins=100)
plt.xlabel("# cells in nhood")
plt.ylabel("# nhoods")
plt.savefig("diff_abundance/nhood_size_histogram.pdf", bbox_inches="tight")
plt.close()
print("Mean neighborhood size:", np.mean(nhood_size))
print("Median neighborhood size:", np.median(nhood_size))

# Count cells in neighbourhoods
# Milo leverages the variation in cell numbers between replicates for the same experimental condition to test for differential abundance.
# We have to count how many cells from each sample are in each neighbourhood.
mmdata = milo.count_nhoods(mmdata, sample_col="sample")

#########################################
# Differential abundance testing with GLM
#########################################
# Create design
experiment_design = "~" + column_name

all_levels = mmdata["rna"].obs[column_name].cat.categories.tolist()
other_levels = [cat for cat in all_levels if cat not in [ref_level, target_level]]
new_order = [ref_level, target_level] + other_levels

# Reorder categories (by default, the last category is taken as the condition of interest)
mmdata["rna"].obs[column_name] = mmdata["rna"].obs[column_name].cat.reorder_categories(new_order)
milo.da_nhoods(mmdata, design=experiment_design, model_contrasts=f"{target_level} - {ref_level}", solver="pydeseq2")
mmdata['milo'].var

# Diagnostic plots on results
with PdfPages(f"diff_abundance/milo_diagnostics.{target_level}_vs_{ref_level}.pdf") as pdf:
    fig = plt.figure(figsize=(10, 5))
    # Panel 1: P-values histogram
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.hist(mmdata["milo"].var.PValue, bins=50)
    ax1.set_xlabel("P-Values")
    # Panel 2: logFC vs -log10(FDR)
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(mmdata["milo"].var.logFC, -np.log10(mmdata["milo"].var.SpatialFDR), ".")
    ax2.set_xlabel("log-Fold Change")
    ax2.set_ylabel("- log10(Spatial FDR)")
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

# Visualize results on embedding
# Building an abstracted graph of neighbourhoods that can be superimposed on the single-cell embedding.
# Each node represents a neighbourhood, and the layout of nodes is determined by the position of the index cell in the UMAP.
# The neighbourhoods displaying singificant DA are colored by their log-Fold Change.
milo.build_nhood_graph(mmdata, basis='X_umap')
plt.figure(figsize=(6, 6))
milo.plot_nhood_graph(mmdata, alpha=0.1, min_size=1, min_logFC=1)
plt.savefig(f"diff_abundance/umap_DA.{target_level}_vs_{ref_level}.pdf", bbox_inches="tight")
plt.close()


# Visualize result by celltype
# Visualize whether DA is particularly evident in certain cell types.
# To do this, each neighbourhood is assigned to a cell type by finding the most abundant cell type within cells in that neighbourhood
milo.annotate_nhoods(mmdata, anno_col="leiden_0.5")
# Mixed celltypes if the fraction of the most abundant celltype is lower than 0.6
mmdata["milo"].var["nhood_annotation"] = mmdata["milo"].var["nhood_annotation"].cat.add_categories("Mixed")
mmdata["milo"].var.loc[mmdata["milo"].var["nhood_annotation_frac"] < 0.6, "nhood_annotation"] = "Mixed"
# Saving table of results
mmdata['milo'].var.to_excel(f"diff_abundance/milo_{target_level}_vs_{ref_level}.xlsx", sheet_name=f"{target_level}_vs_{ref_level}", index=False)

# Diagnostic histogram of neighbourhood 'purity' by celltypes
plt.figure(figsize=(6, 6))
plt.hist(mmdata["milo"].var["nhood_annotation_frac"], bins=30)
plt.xlabel("celltype fraction")
plt.savefig("diff_abundance/hist_nhoods_celltypes_purity.pdf", bbox_inches="tight")
plt.close()

# Visualize the fold changes by cell type annotation.
plt.figure(figsize=(10, 6))
milo.plot_da_beeswarm(mmdata, alpha=0.1)
plt.savefig(f"diff_abundance/milo_logFC_by_celltypes.{target_level}_vs_{ref_level}.pdf", bbox_inches="tight")
plt.close()
