#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import glob
import platform

import anndata
import pandas as pd
import scanpy as sc
import numpy as np


def _parse_to_adata(
    sample: str,
    genes_path: str,
    metadata_path: str,
    mtx_path: str
):

    adata = sc.read_mtx(mtx_path)
    genes = pd.read_csv(genes_path)
    metadata = pd.read_csv(metadata_path)

    adata.var_names = genes['gene_id']
    adata.var = pd.DataFrame(index=genes['gene_id'])
    adata.var['gene_name'] = genes.set_index('gene_id')['gene_name'].reindex(adata.var.index)
    adata.var['gene_name'] = adata.var['gene_name'].astype(str)
    adata.var_names_make_unique()

    adata.obs = metadata.copy()
    cols_to_drop = ["genes", "metadata", "input_type", "mtx"]
    adata.obs = adata.obs.drop(columns=cols_to_drop, errors='ignore')

    if "bc_wells" in adata.obs.columns:
        adata.obs.set_index("bc_wells", inplace=True)



    adata.var["feature_types"] = "Gene Expression"
    adata.var["genome"] = genes.set_index('gene_id')['genome'].reindex(adata.var.index).astype(str)
    adata.var['genome'] = adata.var['genome'].astype(str)



    # reorder columns for 10x mtx files
    adata.var = adata.var[["gene_name","feature_types", "genome"]]
    print(f"Created AnnData for {sample}: {adata.shape[0]} cells × {adata.shape[1]} genes")


    return adata


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.

    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def dump_versions():
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "pandas": pd.__version__,
            "anndata": anndata.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))


def input_to_adata(
    sample: str,
    genes_path: str,
    metadata_path: str,
    mtx_path: str,
    output: str
):
    print(f"Reading data for sample {sample}")
    print(f"  - mtx: {mtx_path}")
    print(f"  - genes: {genes_path}")
    print(f"  - metadata: {metadata_path}")

    # open main data
    adata = _parse_to_adata(
        sample=sample,
        genes_path=genes_path,
        metadata_path=metadata_path,
        mtx_path=mtx_path,
    )

    # write results
    adata.write_h5ad(f"{output}")
    print(f"Wrote h5ad file to {output}")

    # dump versions
    dump_versions()

    return adata



#
# Run main script
#

# create the directory with the sample name
os.makedirs("${meta.id}", exist_ok=True)



# input_type comes from NF module
adata = input_to_adata(
    mtx_path="${inputs}",
    genes_path="${meta.genes}",
    metadata_path="${meta.metadata}",
    output="${meta.id}_${meta.input_type}_matrix.h5ad",
    sample="${meta.id}",
)
