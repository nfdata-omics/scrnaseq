#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================
# Liana works with the log1p-trasformed counts matrix, so log1p transform the data before running the method.

import warnings
import argparse                                             # command line arguments parser
import os                                                   # filesystem utilities
import pathlib                                              # library for handle filesystem paths
import matplotlib
import matplotlib.pyplot as plt                             # plotting library
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc                                         # single-cell data processing
import liana as li                                          # cell-cell interaction inference
from liana.method import logfc, singlecellsignalr,natmi
from liana.method import connectome,cellphonedb,logfc
from liana.method import cellchat,geometric_mean,scseqcomm
import mudata as mu
import sys
import importlib
import importlib.metadata
import yaml


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

warnings.filterwarnings("ignore")
# Set to not interactive
matplotlib.use("Agg")

# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function perform cell-cell interaction analysis using the LIANA framework.
    It takes as input a combined count matrix in h5 format, processes the data, and outputs the results of the cell-cell interaction analysis.
    """

    # --------------------------------------------------------------------------------------------------------------------
    #                                          LIBRARY CONFIG
    # --------------------------------------------------------------------------------------------------------------------

    sc.settings.verbosity = (
        3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    )
    sc.logging.print_header()

    # --------------------------------------------------------------------------------------------------------------------
    #                                          INPUT FROM COMMAND LINE
    # --------------------------------------------------------------------------------------------------------------------

    # Define command line arguments with argparse
    parser = argparse.ArgumentParser(
        description="Cell-cell interaction analysis with LIANA framework",
        epilog="This function performs cell-cell interaction analysis using the LIANA framework. "
            "It takes as input a combined count matrix in h5mu format, processes the data, "
            "and outputs the results of the cell-cell interaction analysis.",
    )
    parser.add_argument(
        "-ad",
        "--input-h5mu-file",
        metavar="H5MU_INPUT_FILES",
        type=pathlib.Path,
        dest="input_h5mu_files",
        required=True,
        help="paths of existing count matrix files in h5 format (including file names)"
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        default="cellphonedb",
        help="method to use for cell-cell interaction analysis (default: cellphonedb)"
    )
    parser.add_argument(
        "-resource",
        "--resource-name",
        type=str,
        default="consensus",
        help="resource to use for cell-cell interaction analysis "
            "(default: consensus for human data, while mouseconsensus is for mouse data)"
    )
    parser.add_argument(
        "-o",
        "--out",
        metavar="H5MU_OUTPUT_FILE",
        type=pathlib.Path,
        default="matrix.ccc.h5mu",
        help="name of the output h5ad file after cell-cell interaction analysis (default: matrix.ccc.h5mu)"
    )
    parser.add_argument(
        "-r",
        "--results",
        type=pathlib.Path,
        default=pathlib.Path("./"),
        help="directory to save the results files (default is the current directory)"
    )
    parser.add_argument(
        "--resolution",
        type=float,
        help="Resolution for the cell-cell interaction analysis"
    )
    parser.add_argument(
        "--versions-dict",
        type=str,
        help="Return dictionary of versions used by the module and exit"
    )

    args = parser.parse_args()

    if args.versions_dict:
        lib_list = ['liana', 'cellphonedb', 'matplotlib', 'scanpy', 'mudata']
        print(versions_yaml(args.versions_dict, lib_list))
        sys.exit(0)


    # --------------------------------------------------------------------------------------------------------------------
    #                                 DEFINE SAMPLES AND MTX PATHS
    # --------------------------------------------------------------------------------------------------------------------

    input_h5mu_file = args.input_h5mu_files
    method = args.method
    resource = args.resource_name
    output_dir = os.path.join(os.getcwd(),f'cell_cell_interaction_final/{method}_{resource}')
    resolution = "leiden_" + str(args.resolution)


    # --------------------------------------------------------------------------------------------------------------------
    #                                 READ H5MU FILES
    # --------------------------------------------------------------------------------------------------------------------

    # Read folders with the MTX combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # Read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder ... ", end="")
    mudata = mu.read_h5mu(input_h5mu_file)
    adata = mudata.mod['normalized_gex']
    print("Done!")
    print(f"Count matrix for combined samples has {adata.shape[0]} cells and {adata.shape[1]} genes/ab")

    # --------------------------------------------------------------------------------------------------------------------
    #                                 PRE-PROCESSING
    # --------------------------------------------------------------------------------------------------------------------

    sc.pl.umap(adata, color=[resolution], frameon=False)

    # Convert the data to memory-mapped format for efficient processing
    adata = adata.to_memory()

    # --------------------------------------------------------------------------------------------------------------------
    #                                 CELL-CELL INTERACTION ANALYSIS WITH LIANA USING CELLPHONEDB METHOD
    # --------------------------------------------------------------------------------------------------------------------

    print(
        "\n===== CELL-CELL INTERACTION ANALYSIS WITH LIANA USING CHOSEN METHOD  ====="
    )

    methods_dict = {
        "singlecellsignalr": singlecellsignalr,
        "connectome": connectome,
        "cellphonedb": cellphonedb,
        "natmi": natmi,
        "logfc": logfc,
        "cellchat": cellchat,
        "geometric_mean": geometric_mean,
        "scseqcomm": scseqcomm
    }

    method_name = args.method

    if method_name not in methods_dict:
        raise ValueError(
            f"Method '{method_name}' is not valid. Choose from: {list(methods_dict.keys())}"
    )
    method = methods_dict[method_name]

    # Run the chosen method, extracting the top interactions based on the specified resource and parameters.
    method(
        adata,
        groupby=resolution,
        resource_name=resource,    # by default the resource uses HUMAN gene symbols, if you are working with mouse data, you can specify the resource_name as 'mouseconsensus'
        expr_prop=0.1,             # Minimum expression proportion for the ligands and receptors (+ their subunits) in the corresponding cell identities
        return_all_lrs=False,      # Ligand-receptor pairs that pass the expr_prop threshold
        use_raw=False,             # whether to use the log1p-transformed data for the analysis, by default it uses the log1p-transformed data, which is more suitable for the analysis.
        key_added = f"{method_name}_res",
        verbose=True,
        seed=42,
    )

    # By default, any interactions in which either entity is not expressed in above 10% of cells per cell type is considered as a false positive,
    # under the assumption that since CCC occurs between cell types, a sufficient proportion of cells within should express the genes.
    print(adata.uns[f"{method_name}_res"].head())
    res = adata.uns[f"{method_name}_res"]
    res.to_excel(os.path.join(output_dir, f"{method_name}_results.xlsx"),index=False)

    # --------------------------------------------------------------------------------------------------------------------
    #                           VISUALIZE PLOTS OF SIGNIFICANT INTERACTIONS FOR CHOSEN METHOD
    # --------------------------------------------------------------------------------------------------------------------

    # Define plotting config for each method
    plot_config = {
        "cellphonedb": {
            "orderby": "lr_means",
            "orderby_ascending": True,
            "filter_fun": lambda x: x["cellphone_pvals"] <= 0.05,
            "size": "cellphone_pvals",
            "colour": "lr_means",
            "fill": "means",
            "label": "props",
        },
        "connectome": {
            "orderby": "expr_prod",
            "orderby_ascending": False,
            "filter_fun": lambda x: x["scaled_weight"] <= 0.05,
            "size": "scaled_weight",
            "colour": "expr_prod"
        },
        "cellchat": {
            "orderby": "lr_probs",
            "orderby_ascending": False,
            "filter_fun": lambda x: x["cellchat_pvals"] <= 0.05,
            "size": "cellchat_pvals",
            "colour": "lr_probs"
        }
    }

    cfg = plot_config[method_name]

    # Create the dotplot with filtered interactions
    print("\nVisualizing dotplot of significant interactions for chosen method")
    # Define figure size (larger for readability)
    plot_name = f"dotplot_cell-cell_interaction_{method_name}_{resolution}.pdf"
    fig1 = li.pl.dotplot(
        adata=adata,
        colour=cfg["colour"],
        size=cfg["size"],
        inverse_size=True,                                           # small p-values have large sizes
        orderby=cfg["orderby"],
        orderby_ascending=True,
        top_n=20,
        uns_key=f"{method_name}_res",                                # use dynamic key_added
        filter_fun=cfg["filter_fun"],
        size_range=(1, 6),
        return_fig=True,
    )

    # Save the figure
    fig1.save(
        os.path.join(output_dir, plot_name),
        width=60,
        height=20,
        limitsize=False
    )



    # Create the tileplot with filtered interactions statistics of ligands and receptors across the interacting cell types, which can be used to identify the cell types that are most likely to be involved in the interactions.
    plot_name = f"tileplot_cell-cell_interaction_{method_name}_{resolution}.pdf"
    fig2 = li.pl.tileplot(
        adata=adata,
        fill=cfg["fill"],
        label=cfg["label"],
        label_fun=lambda x: f"{x:.2f}",
        top_n=20,
        orderby=cfg["orderby"],
        orderby_ascending=True,
        uns_key=f"{method_name}_res",  # uns_key to use, default is 'liana_res'
        source_title="Ligand",
        target_title="Receptor",
        filter_fun=cfg["filter_fun"],
        return_fig=True,
    )
    fig2.save(
        os.path.join(output_dir, plot_name),
        width=18,
        height=14,
    )

if __name__ == "__main__":
    main()
