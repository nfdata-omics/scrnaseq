#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."
os.environ["MPLCONFIGDIR"] = "/tmp"

import platform
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

def read_samplesheet(samplesheet):
    df = pd.read_csv(samplesheet)
    df.set_index("sample")
    if 'feature_type' in df.columns:
        df['feature_type'] = df['feature_type'].fillna('unknown').astype(str)
    elif 'sample_type' in df.columns:
        df['sample_type'] = df['sample_type'].fillna('unknown').astype(str)
    else:
        print("Warning: Neither 'feature_type' nor 'sample_type' found.")


    # samplesheet may contain replicates, when it has,
    # group information from replicates and collapse with commas
    # only keep unique values using set()
    df = df.groupby(["sample"]).agg(lambda column: ",".join(set(column.astype(str))))

    return df


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
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))


if __name__ == "__main__":
    # Open samplesheet as dataframe
    df_samplesheet = read_samplesheet("${samplesheet}")

    # find all h5ad and append to dict
    dict_of_h5ad = {}

    for path in Path(".").rglob("*.h5ad"):
        adata_tmp = sc.read_h5ad(path)

        if "feature_types" in adata_tmp.var.columns:
            adata_tmp = adata_tmp[:, adata_tmp.var["feature_types"] != "Peaks"].copy()

        key = str(path).replace("_matrix.h5ad", "")
        dict_of_h5ad[key] = adata_tmp

    # concat h5ad files
    adata = ad.concat(dict_of_h5ad, label="sample",join="outer", index_unique="_")

    # grab all var DataFrames from dictionary
    all_var = [x.var for x in dict_of_h5ad.values()]
    # concatenate them
    all_var = pd.concat(all_var, join="outer")
    # remove duplicates
    all_var = all_var[~all_var.index.duplicated()]
    adata.var = all_var.loc[adata.var_names]

    # merge with data.frame, on sample information
    adata.obs['sample'] = adata.obs['sample'].str.replace('_filtered', '', regex=False).str.strip()
    adata.obs = adata.obs.join(df_samplesheet, on="sample", how="left").astype(str)
    print(adata.obs)
    adata.write_h5ad("${meta.id}_${meta.input_type}_matrix.h5ad")

    print("Wrote h5ad file to {}".format("${meta.id}_${meta.input_type}_matrix.h5ad"))

    # dump versions
    dump_versions()
