#!/usr/bin/env python

###################################################################
###################### Pseudobulk analysis ########################
###################################################################

# Importing libraries
import warnings
import argparse                     # command line arguments parser
import os                           # filesystem utilities
import pathlib                      # library for handling filesystem paths
import matplotlib.pyplot as plt     # library for visualization
import pandas as pd                 # library for data analysis and manipulation
import scanpy as sc                 # single-cell data processing
import decoupler as dc
import pertpy as pt
import mudata as md
import sys
import importlib
import importlib.metadata
import yaml

warnings.filterwarnings("ignore")

def versions_yaml(process_name, list_of_libs=None):
    """
    Generate YAML formatted string with versions of relevant libraries.

    Parameters
    ----------
    process_name : str
        Process name to use as key in the versions dictionary.
    list_of_libs : list of str, optional
        List of specific library names to include in the versions dictionary.

    Returns
    -------
    str
        YAML formatted string containing library versions and Python version.
    """

    versions = {}
    versions[process_name] = {}

    versions[process_name]['python'] = f"{sys.version_info.major}" \
        f".{sys.version_info.minor}.{sys.version_info.micro}"

    for lib in list_of_libs:
        try:
            version = importlib.metadata.version(lib)
        except importlib.metadata.PackageNotFoundError:
            try:
                module = importlib.import_module(lib)
                version = getattr(module, '__version__', 'unknown')
            except (ImportError, AttributeError):
                version = None
        if version is not None:
            versions[process_name][lib] = version

    return yaml.dump(versions)

#########################################
# Define and parse CLI arguments
#########################################
parser = argparse.ArgumentParser(description='Pseudobulk')

parser.add_argument("--mdata", type=str,
                    help="Path to input MuData object (h5mu file)")

parser.add_argument("--group_column", type=str,
                    help="Column in the metadata to group by for pseudobulk aggregation (e.g., genotype, condition)")

parser.add_argument("--resolution", type=float,
                    help="Clustering resolution to use for pseudobulk aggregation")

parser.add_argument("--comparisons", type=str,
                    help="Comparisons to be performed following column:target:reference (e.g. group:treated:control)")

parser.add_argument("--versions-dict", type=str,
                    help="Return dictionary of versions used by the module and exit")

args = parser.parse_args()

if args.versions_dict:
   lib_list = ['pandas', 'pertpy', 'mudata', 'argparse', 'matplotlib', 'scanpy', 'decoupler']
   print(versions_yaml(args.versions_dict, lib_list ))
   sys.exit(0)

#########################################
# Get column_name, target_level and
# ref_level from comparisons
#########################################
comparisons = args.comparisons.split(':')

column_name = comparisons[0]
target_level = comparisons[1]
ref_level = comparisons[2]

# Create output directory
os.makedirs('pseudobulk', exist_ok=True)

# Import object
mdata = md.read(args.mdata)
gex = mdata.mod['gex']

# Handle gene names (parse version)
gex.var['gene_symbols'] = gex.var['gene_symbols'].astype('str')
gex.var_names = gex.var["gene_symbols"]
gex.var_names_make_unique()

# Check if the selected resolution is present in object metadata
set_res = args.resolution
cluster_key = "leiden_{}".format(set_res)
if cluster_key not in gex.obs.keys():
    print(f"Clustering at resolution {set_res} not found. Exiting.")
    sys.exit(1)

# Extract pseudobulk profiles within each group and sample of interest
ps = pt.tl.PseudobulkSpace()

# Set variable with column name to perform pseudobulking on
group_column = args.group_column

gex.obs[f'sample_{set_res}'] = (gex.obs['sample'].astype(str) + "_" + gex.obs[cluster_key].astype(str))

# Calculate pseudobulk
pdata = ps.compute(gex,
                   target_col=f"sample_{set_res}",
                   groups_col=group_column,
                   layer_key="count",
                   mode="sum",
                   min_cells=10,
                   min_counts=1000)

# Plotting counts vs ncells
fig, ax = plt.subplots(figsize=(40,20))
ps.plot_psbulk_samples(pdata,
                       groupby=[f"leiden_{set_res}", group_column],
                       figsize=(12, 6),
                       return_fig=True)
plt.savefig(os.path.join(os.getcwd(),f"pseudobulk/pseudobulk_counts_{set_res}.pdf"),
            bbox_inches='tight',
            dpi=300)
plt.close()

# PCA
pdata.layers["counts"] = pdata.X.copy()
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.pp.pca(pdata)
dc.swap_layer(pdata, "counts", inplace=True)

# Plotting
variables = ["sample", f"leiden_{set_res}", group_column]
output_dir = os.path.join(os.getcwd(), 'pseudobulk')
os.makedirs(output_dir, exist_ok=True)

for var in variables:
    fig, ax = plt.subplots(figsize=(40, 20))
    sc.pl.pca(pdata, color=var, size=300, show=False)
    file_path = os.path.join(output_dir, f'PCA_{var}.pdf')
    plt.savefig(file_path, bbox_inches='tight')
    plt.close()

# PCA variance ratio
fig, ax = plt.subplots(figsize=(40, 20))
sc.pl.pca_variance_ratio(pdata, show=False)
file_path = os.path.join(output_dir, 'PCA_variance_ratio.pdf')
plt.savefig(file_path, bbox_inches='tight')
plt.close()

# Get clusters at selected resolution
clusters = gex.obs[f"leiden_{set_res}"].astype(str).unique().tolist()

for cl in clusters:
    print(f"\nProcessing cluster {cl}")

    # Subset cells to one cluster
    gex_cl = gex[gex.obs[f"leiden_{set_res}"].astype(str) == cl].copy()

    # Skip small clusters
    if gex_cl.n_obs < 50:
        print(f"Skipping cluster {cl} (less than 50 cells)")
        continue

    # Calculate pseudobulk for cluster
    pdata = ps.compute(
        gex_cl,
        target_col="sample",
        groups_col=group_column,
        layer_key="count",
        mode="sum",
        min_cells=10,
        min_counts=1000
    )

    pdata.layers["counts"] = pdata.X.copy()

    contrast_groups = [ref_level, target_level]
    groups = pdata.obs[group_column].unique()
    if not set(contrast_groups).issubset(groups):
        print(f"Skipping cluster {cl} (missing groups)")
        continue

    counts = pdata.obs.loc[
        pdata.obs[group_column].isin(contrast_groups),
        group_column
    ].value_counts()

    if (counts < 2).any():
        print(f"Skipping cluster {cl} (insufficient replicates per group)")
        continue

    # Export data for analysis with DESeq2
    export_dir = os.path.join(os.getcwd(), "pseudobulk", "export_deseq2", cl)
    os.makedirs(export_dir, exist_ok=True)

    counts_df = pd.DataFrame(
        pdata.layers["counts"].T,
        index=pdata.var_names,
        columns=pdata.obs_names
    )
    coldata_df = pdata.obs.copy()
    coldata_df.index.name = "sample"

    counts_df.to_csv(os.path.join(export_dir, f"counts_cl_{cl}.tsv"), sep="\t")
    coldata_df.to_csv(os.path.join(export_dir, f"coldata_cl_{cl}.tsv"), sep="\t")

    print(f"Saved cluster {cl}")
