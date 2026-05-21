#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import pathlib                      # library for handle filesystem paths
from pathlib import Path
import glob
import scanpy as sc                 # single-cell data processing
import scirpy as ir                 # single-cell AIRR-data
import anndata as ad                # store annotated matrix as anndata object
import pandas as pd                 # data analysis and manipulation tool
import os                           # misc operating system interfaces
import mudata as md
import matplotlib.pyplot as plt     # plotting library
import numpy as np                  # scientific computing with Python
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
    This function concatenates csv files from vdj modality.
    """

# --------------------------------------------------------------------------------------------------------------------
#                                          LIBRARY CONFIG
# --------------------------------------------------------------------------------------------------------------------

    sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

#Define command line arguments with argparse
    parser = argparse.ArgumentParser(prog='QC_filter_vdj', usage='%(prog)s [options]', description = "QC metrics and filtering for vdj data",
                                    epilog = "This function calculates common quality control (QC) metrics for each sample for vdj modality, inspects QC plots for each sample.",
                                    )
    parser.add_argument('-ad','--input-h5mu-combined',metavar= 'H5MU_INPUT_FILES', type=pathlib.Path, dest='input_h5mu_files',
                        required=True, help="paths of existing matrix files in h5mu format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5MU_OUTPUT_FILE', type=pathlib.Path, default="matrix.qc_vdj.h5mu",
                        help="path and name of the output h5mu file")
    parser.add_argument('-csv', '--csv_out', metavar='QUALITY_CONTROL', default="group_abundance.csv",
                        help="path and name of csv table with ranked marker genes for each cluster and resolution")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'),
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5MU FILES =====")
    input_h5mu_file = args.input_h5mu_files
    output_csv= Path(args.csv_out)
    output =args.out


    # print info on the available matrices
    print("Reading combined matrix from the following file:")
    print(f"-File {input_h5mu_file}")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5MU FILES
# --------------------------------------------------------------------------------------------------------------------

    # Read folders with the combined count matrice and store datasets in a dictionary
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print(f"\nProcessing MuData object in folder {input_h5mu_file} ... ", end ='')

    mdata= md.read(input_h5mu_file)
    print("Done!")
    print(f"MuData matrix for combined samples has {mdata.shape[0]} cells and {mdata.shape[1]} genes/ab")
    print(mdata)
# --------------------------------------------------------------------------------------------------------------------
#                                 VDJ MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== VDJ MODALITY DATA =====")

    if "airr" not in mdata.mod:
        raise ValueError("No 'airr' modality found in MuData")

    vdj = mdata.mod["airr"]
    print(f"VDJ data found with {vdj.n_obs} cells and {vdj.n_vars} features.")
    print(vdj.obsm)

# --------------------------------------------------------------------------------------------------------------------
#                                 GEX MODALITY DATA
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== GEX MODALITY DATA =====")
    # Check if 'GEX' exists in mdata.mod
    if "gex" not in mdata.mod:
        raise ValueError("No 'gex' modality found in MuData")

    gex = mdata.mod["gex"]
    print(f"GEX data found with {gex.n_obs} cells and {gex.n_vars} features.")
    print(gex.obs)


# --------------------------------------------------------------------------------------------------------------------
#                           CREATING CHAIN INDICES
# --------------------------------------------------------------------------------------------------------------------
    print("\n===== CREATING CHAIN INDICES =====")
    n_before = sum(len(x) for x in vdj.obsm["airr"])
    print(f"Number of VDJ entries before creating chain indices: {n_before}")
    # Create chain indices
    ir.pp.index_chains(vdj,filter = ["productive","require_junction_aa"])
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           TCR QUALITY CONTROL
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== TCR QUALITY CONTROL =====")
    # Compute basic quality control metrics
    ir.tl.chain_qc(vdj)
    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           GROUP ABUNDACE
# --------------------------------------------------------------------------------------------------------------------

    all_abundance = []

    print("\n===== GROUP ABUNDANCE =====")

    for sample in vdj.obs["sample"].unique():

        print(f"\nProcessing sample {sample}")

        ad = vdj[vdj.obs["sample"] == sample]

        # receptor_type
        df_type = ir.tl.group_abundance(
            ad,
            groupby="receptor_type",
            target_col="sample",
            fraction=False
        )

        df_type = (
            df_type
            .stack()
            .reset_index()
        )

        df_type.columns = ["category", "sample", "count"]
        df_type["metric"] = "receptor_type"

        all_abundance.append(df_type)

        # receptor_subtype
        df_sub = ir.tl.group_abundance(
            ad,
            groupby="receptor_subtype",
            target_col="sample",
            fraction=False
        )

        df_sub = (
            df_sub
            .stack()
            .reset_index()
        )

        df_sub.columns = ["category", "sample", "count"]
        df_sub["metric"] = "receptor_subtype"

        all_abundance.append(df_sub)

        # chain_pairing
        df_chain = ir.tl.group_abundance(
            ad,
            groupby="chain_pairing",
            target_col="sample",
            fraction=False
        )

        df_chain = (
            df_chain
            .stack()
            .reset_index()
        )

        df_chain.columns = ["category", "sample", "count"]
        df_chain["metric"] = "chain_pairing"

        all_abundance.append(df_chain)


    combined_df = pd.concat(all_abundance, ignore_index=True)

    # reorder columns
    combined_df = combined_df[
        ["sample", "metric", "category", "count"]
        ]

    output_file = output_csv.parent / "VDJ_abundance_all_samples.csv"

    combined_df.to_csv(output_file, index=False)

    print(f"\nSaved combined abundance table at {output_file}")


# --------------------------------------------------------------------------------------------------------------------
#                           STATISTICS ON THE NUMBER OF CHAINS PER CELL
# --------------------------------------------------------------------------------------------------------------------

    types = ["single pair","extra VJ", "extra VDJ", "two full chains"]


    rows = []

    for sample in vdj.obs["sample"].unique():
        n_total = (vdj.obs["sample"] == sample).sum()
        n = np.sum(
            (vdj.obs["sample"] == sample)
            & (vdj.obs["chain_pairing"].isin(types))
        )
        fraction = n / n_total if n_total > 0 else 0

        rows.append({
            "sample": sample,
            "n_total": n_total,
            "n_extra": n,
            "fraction": fraction
        })

    output_file = output_csv.parent / "chain_pairing_stats_all_samples.csv"
    pd.DataFrame(rows).to_csv(output_file, index=False)
    print(f"Saved chain pairing statistics for all samples at {output_file}")


# --------------------------------------------------------------------------------------------------------------------
#                           GROUP ABUNDACE PLOTTING
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== GROUP ABUNDACE PLOTTING =====")
    pdf_type = os.path.join(args.results, "vdj_receptor_type.pdf")
    pdf_subtype = os.path.join(args.results, "vdj_receptor_subtype.pdf")
    pdf_chain = os.path.join(args.results, "vdj_chain_pairing.pdf")

    with PdfPages(pdf_type) as pdf_type_file, PdfPages(pdf_subtype) as pdf_subtype_file, PdfPages(pdf_chain) as pdf_chain_file:
        for sample in vdj.obs["sample"].unique():

            print(f"\nProcessing sample {sample}")
            sample_data = vdj[vdj.obs["sample"] == sample]

            # Receptor Type
            ir.pl.group_abundance(
                sample_data,
                groupby="receptor_type",
                target_col="sample")
            fig = plt.gcf()
            fig.set_size_inches(10, 8)
            fig.suptitle(f"Receptor Type - Sample: {sample}", fontsize=16, fontweight='bold', y=1.02)
            pdf_type_file.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            print(f"Added receptor type page for sample {sample}")

            # Receptor Subtype
            ir.pl.group_abundance(
                sample_data,
                groupby="receptor_subtype",
                target_col="sample")
            fig = plt.gcf()
            fig.set_size_inches(10, 8)
            fig.suptitle(f"Receptor Subtype - Sample: {sample}", fontsize=16, fontweight='bold', y=1.02)
            pdf_subtype_file.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            print(f"Added receptor subtype page for sample {sample}")

            # Chain Pairing
            ir.pl.group_abundance(
                sample_data,
                groupby="chain_pairing",
                target_col="sample")
            fig = plt.gcf()
            fig.set_size_inches(10, 8)
            fig.suptitle(f"Chain Pairing - Sample: {sample}", fontsize=16, fontweight='bold', y=1.02)
            pdf_chain_file.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            print(f"Added chain pairing page for sample {sample}")

    print(f"Saved receptor type PDF at {pdf_type}")
    print(f"Saved receptor subtype PDF at {pdf_subtype}")
    print(f"Saved chain pairing PDF at {pdf_chain}")



# --------------------------------------------------------------------------------------------------------------------
#                          MATCH VDJ METRICS TO RNA MODALITY
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== MATCH VDJ METRICS TO RNA MODALITY =====")
    # Match VDJ metrics to RNA modality based on cells with productive receptors

    # Clean sample names and barcodes
    vdj.obs["sample_clean"] = (
        vdj.obs["sample"].astype(str)
        .str.replace("_cellbender_filter", "", regex=False)
        .str.replace("_filtered", "", regex=False)
        .str.replace("_parse", "", regex=False)
        .str.strip()
    )

    gex.obs["sample_clean"] = (
        gex.obs["sample"].astype(str)
        .str.replace("_cellbender_filter", "", regex=False)
        .str.replace("_filtered", "", regex=False)
        .str.replace("_parse", "", regex=False)
        .str.strip()
    )


    vdj.obs["barcode_clean"] = (
        vdj.obs_names.astype(str)
        .str.split("_")
        .str[0]
    )

    gex.obs["barcode_clean"] = (
        gex.obs_names.astype(str)
        .str.split("_")
        .str[0]
    )

    vdj.obs["match_key"] = (
        vdj.obs["sample_clean"].astype(str)
        + "_"
        + vdj.obs["barcode_clean"].astype(str)
    )
    gex.obs["match_key"] = (
        gex.obs["sample_clean"].astype(str)
        + "_"
        + gex.obs["barcode_clean"].astype(str)
    )

    match_stats = []

    types = ["single pair", "extra VJ", "extra VDJ", "two full chains"]

    vdj_samples = set(vdj.obs["sample_clean"].unique())
    gex_samples = set(gex.obs["sample_clean"].unique())
    common_samples = vdj_samples & gex_samples

    print("Samples in VDJ:", sorted(vdj_samples))
    print("Samples in GEX:", sorted(gex_samples))
    print("Common samples:", sorted(common_samples))
    print(f"VDJ sample count: {len(vdj_samples)}")
    print(f"GEX sample count: {len(gex_samples)}")
    print(f"Common sample count: {len(common_samples)}")

    all_vdj_keys = vdj.obs.loc[
        vdj.obs["chain_pairing"].isin(types),
        "match_key"
    ].astype(str)
    all_gex_keys = gex.obs.loc[:, "match_key"].astype(str)
    common_keys = np.intersect1d(all_vdj_keys, all_gex_keys)
    print(f"Common sample+barcode keys between VDJ and GEX: {len(common_keys)}")

    for sample in sorted(common_samples):
        vdj_cells = vdj.obs.loc[
            (vdj.obs["sample_clean"] == sample) &
            (vdj.obs["chain_pairing"].isin(types)),
            "match_key"
        ].astype(str)

        gex_cells = gex.obs.loc[
            gex.obs["sample_clean"] == sample,
            "match_key"
        ].astype(str)

        n_vdj = len(vdj_cells)
        n_gex = len(gex_cells)

        matching_cells = np.intersect1d(vdj_cells, gex_cells)
        n_match = len(matching_cells)

        match_stats.append({
            "sample": sample,
            "vdj_cells": n_vdj,
            "gex_cells": n_gex,
            "matched_with_gex": n_match,
            "fraction_matched_on_total_vdj_cells": n_match / n_vdj if n_vdj > 0 else 0
        })

        print(f"{sample}: {n_match}/{n_vdj} VDJ cells match GEX ({n_match/n_vdj:.2%})")

    match_df = pd.DataFrame(match_stats)
    output_file = output_csv.parent / "vdj_gex_match_all_samples.csv"
    match_df.to_csv(output_file, index=False)

    print(f"\nSaved VDJ-GEX match statistics at {output_file}")


# --------------------------------------------------------------------------------------------------------------------
#                           GROUP ABUNDACE ON INFERRED DONORS
# --------------------------------------------------------------------------------------------------------------------

    # productive_types = ["single pair", "extra VJ", "extra VDJ", "two full chains", "multichain"]
    # vdj_productive = vdj[vdj.obs["chain_pairing"].isin(productive_types)]

    # vdj_matched = vdj_productive[vdj_productive.obs_names.isin(gex.obs["barcode_clean"])]

    # donor_map = gex.obs[["barcode_clean", "Inferred_donor"]].set_index("barcode_clean")
    # vdj_matched.obs["Inferred_donor"] = vdj_matched.obs_names.map(donor_map["Inferred_donor"])


    # for donor in vdj_matched.obs["Inferred_donor"].dropna().unique():
    #     df = ir.tl.group_abundance(
    #         vdj_matched[vdj_matched.obs["Inferred_donor"] == donor],
    #         groupby="chain_pairing",
    #         target_col="Inferred_donor",
    #         fraction=False
    #     )
    #     output_file = output_csv.parent / f"ChainPairing_donor-{donor}.csv"
    #     df.to_csv(output_file, index=True)
    #     print(f"Saved table for donor {donor} at {output_file}")

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------

    mdata.mod["airr"] = vdj
    #mdata.mod["gex"] = gex

    #for mod in ["gex", "airr"]:
    #    for col in ["barcode_clean", "sample_clean"]:
    #        if col in mdata.mod[mod].obs.columns:
    #            mdata.mod[mod].obs.drop(columns=[col], inplace=True)


    print("\n===== SAVING OUTPUT FILE =====")
    print(f"Saving h5mu data to file {output}")
    mdata.write(output)
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
