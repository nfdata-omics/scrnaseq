#!/usr/bin/env python

# Set numba cache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import platform
from pathlib import Path

import anndata
import pandas as pd
import scanpy as sc


def read_10x_matrix(
    input_data: str,
    sample: str,
):
    input_path = Path(input_data)

    if input_path.is_file() and input_path.suffix.lower() == ".h5":
        adata = sc.read_10x_h5(input_path, gex_only=False)
    elif input_path.is_dir():
        adata = sc.read_10x_mtx(input_path, gex_only=False)
    else:
        raise ValueError(
            f"Unsupported matrix input: {input_path}. "
            "Expected an H5 file or a complete MEX directory."
        )

    if "gene_ids" not in adata.var:
        raise ValueError(f"Matrix input {input_path} does not contain gene identifiers.")

    adata.var["gene_symbols"] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)
    adata.obs["sample"] = sample

    # Keep the same leading columns across H5 and MEX while preserving optional
    # metadata, such as the genome column available in Cell Ranger H5 files.
    preferred_columns = ["gene_symbols", "feature_types", "genome"]
    ordered_columns = [
        column for column in preferred_columns if column in adata.var.columns
    ]
    ordered_columns += [
        column for column in adata.var.columns if column not in ordered_columns
    ]
    adata.var = adata.var[ordered_columns]

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
    input_data: str,
    output: str,
    sample: str,
):
    print(f"Reading in {input_data}")

    adata = read_10x_matrix(input_data, sample)

    # standard format
    # index are gene IDs and symbols are a column
    adata.var["gene_versions"] = adata.var.index
    adata.var.index = adata.var["gene_versions"].str.split(".").str[0].values
    adata.var_names_make_unique()

    # write results
    adata.write_h5ad(output)
    print(f"Wrote h5ad file to {output}")

    # dump versions
    dump_versions()


#
# Run main script
#

input_to_adata(
    input_data="${matrix_input}",
    output="${meta.id}_${meta.input_type}_matrix.h5ad",
    sample="${meta.id}",
)
