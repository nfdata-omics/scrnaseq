#!/usr/bin/env python3
import pandas as pd
import anndata as ad
import mudata as md
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import mllmcelltype as mct
from datetime import datetime
import time
import os
from dotenv import load_dotenv
import argparse
from pathlib import Path

# set script version number
VERSION = "0.0.1"

def main():
    """
    This function performs cell type annotation using mLLMCelltype for multiple Leiden resolutions using various models from OpenRouter, and saves the annotated clusters and UMAP plots.
    """

    # Define command line arguments with argparse
    parser = argparse.ArgumentParser(description='mLLMCelltype')
    parser.add_argument('--input', type=Path, help='Input file (.h5ad or .h5mu) with precomputed Leiden clusters')
    parser.add_argument('--species', type=str, help='Species for annotation')
    parser.add_argument('--tissue', type=str, help='Tissue type for annotation')
    parser.add_argument('--n_genes', type=int, default=10, help='Number of marker genes to use for annotation')
    parser.add_argument('--pval', type=float, default=0.01, help='P-value threshold for marker gene selection')
    parser.add_argument('--lfc', type=float, default=1.5, help='Log fold change threshold for marker gene selection')
    parser.add_argument('--resolutions', type=str, nargs='+', help='Leiden resolutions to do cell type annotation on (e.g. leiden_0.3)')
    parser.add_argument('--umap_embedding', type=str, default='X_umap', help='Name of UMAP embedding in .obsm (default: X_umap)')
    parser.add_argument('--out', metavar='H5MU_OUTPUT_FILE', type=Path, default="matrix.annotated.h5mu", help="name of the output h5mu file after cell annotation")
    args = parser.parse_args()

    # API Keys - FIGURE OUT WHERE TO PUT THIS IN PIPELINE
    load_dotenv()
    openrouter_api_key = os.environ["OPENROUTER_API_KEY"]
    gemini_api_key = os.environ["GEMINI_API_KEY"]
    api_keys = {
        "openrouter": openrouter_api_key,
        "gemini": gemini_api_key,
    }

    # Force mllmcelltype to use OpenRouter for all OpenAI-compatible models
    os.environ['OPENAI_API_BASE'] = 'https://openrouter.ai/api/v1'
    os.environ['OPENAI_BASE_URL'] = 'https://openrouter.ai/api/v1'

    # Model configuration - all models with their providers
    models_config = {
        'openai/gpt-4.1-nano': 'openrouter'
        }

    # Setup
    mct.setup_logging(log_level='ERROR')
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    sc.settings.figdir = "." # Set scanpy figure directory

    # --------------------------------------------------------------------------------------------------------------------
    #                                 READ H5MU FILES
    # --------------------------------------------------------------------------------------------------------------------

    # Load input (auto-detect based on extension)
    print("\n===== READING MUDATA =====")
    if args.input.suffix == '.h5mu':
        mu = md.read_h5mu(args.input)
        adata = mu.mod['gex']
        print(f"Loaded MuData with modalities: {list(mu.mod.keys())}")
    elif args.input.suffix == '.h5ad':
        adata = ad.read_h5ad(args.input)
        print("Loaded AnnData")
    else:
        raise ValueError(f"Unsupported file format: {args.input.suffix}. Use .h5ad or .h5mu")

    print("Available obsm columns:\n", adata.obsm.keys())

    # -------------------------------------------------------------------------------------------------------------------
    #                                 MLLMCELLTYPE ANNOTATION
    # -------------------------------------------------------------------------------------------------------------------

    all_results = []

    for res in args.resolutions:
        print(f"\n=== Processing Leiden resolution:: {res} ===")

        # Calculate marker genes (once per resolution)
        sc.tl.rank_genes_groups(adata, res, method='wilcoxon')
        marker_genes = {}
        for cluster in adata.obs[res].unique():
            genes = sc.get.rank_genes_groups_df(adata, group=cluster, pval_cutoff=args.pval, log2fc_min=args.lfc)['names'].tolist()[:args.n_genes]
            marker_genes[cluster] = genes
        print(f"Found {len(marker_genes)} clusters\n")
        print(f"Marker Genes per Cluster: {marker_genes}\n")

        for model, provider in models_config.items():
            time.sleep(5)
            print(f"Testing model {model} with provider {provider}:")
            model_short = model.split('/')[1].split(':')[0] if '/' in model else model
            try:
                # Retry logic with exponential backoff
                attempts = 1
                annotations = None
                for attempt in range(attempts):
                    try:
                        annotations = mct.annotate_clusters(
                            marker_genes=marker_genes,
                            species=args.species,
                            tissue=args.tissue,
                            provider=provider,
                            model=model,
                            api_key=api_keys[provider]
                        )
                        break  # Exit loop if successful
                    except Exception as e:
                        if attempt < attempts - 1:
                            wait_time = 2 ** attempt # exponential backoff
                            print(f"ERROR: {str(e)}. Retrying in {wait_time} seconds...")
                            time.sleep(wait_time)
                        else:
                            print(f"ERROR: {str(e)}. No more attempts left.")
                            raise

                print("Annotations:")
                for cluster, cell_type in sorted(annotations.items()):
                    print(f"  Cluster {cluster}: {cell_type}")
                print()

                # Store cluster annotations
                for cluster, celltype in annotations.items():
                    all_results.append({
                        'resolution': res,
                        'model': model_short,
                        'cluster': cluster,
                        'cell_type': celltype
                    })

                # --------------------------------------------------------------------------------------------------------------------
                #                           ADD ANNOTATIONS BACK TO ANNDATA
                # --------------------------------------------------------------------------------------------------------------------
                col_name = f"celltype_{res}_{model_short}"
                col_name_with_cluster = f"{col_name}_with_cluster"

                # Map cluster IDs to cell types
                adata.obs[col_name] = adata.obs[res].astype(str).map(annotations)

                # Create combined label with cluster number and cell type
                adata.obs[col_name_with_cluster] = (
                    adata.obs[res].astype(str) + ": " + adata.obs[col_name].astype(str)
                )

                print(f"Annotating with {model}...")
                print(f"  Unique cell types: {adata.obs[col_name].nunique()}")

                # --------------------------------------------------------------------------------------------------------------------
                #                           PLOT UMAP WITH ANNOTATIONS
                # --------------------------------------------------------------------------------------------------------------------
                plot_name = f'umap_{res}_{model_short}_{timestamp}'
                try:
                    sc.pl.embedding(
                        adata,
                        basis=args.umap_embedding,
                        color=col_name_with_cluster,
                        title=f'{model} - {res}',
                        show=False
                    )

                    plt.savefig(os.path.join("./", f"{plot_name}.pdf"), bbox_inches='tight', dpi=300)
                    plt.close()
                    print(f"Saved UMAP plot as: {os.path.join("./", f'{plot_name}.pdf')}\n")
                except Exception as plot_error:
                    print(f"Warning: Could not save UMAP plot: {plot_error}\n")
            except Exception as e:
                print(f"ERROR: {str(e)}\n")

    # --------------------------------------------------------------------------------------------------------------------
    #                           SAVE ANNOTATIONS TO CSV
    # --------------------------------------------------------------------------------------------------------------------
    results_df = pd.DataFrame(all_results)
    results_file = f'./annotated_clusters_{timestamp}.csv'
    results_df.to_csv(results_file, index=False)
    print(f"Saved annotation results to {results_file}")

    # --------------------------------------------------------------------------------------------------------------------
    #                           SAVE INPUT PARAMETERS TO TXT
    # --------------------------------------------------------------------------------------------------------------------
    with open(f'./parameters_{timestamp}.txt', 'w') as f:
        f.write(f"Models: {list(models_config.keys())}\np-value cutoff: {args.pval}\nlog2fc minimum: {args.lfc}\nNumber of marker genes: {args.n_genes}\nResolutions: {args.resolutions}\nSpecies: {args.species}\nTissue: {args.tissue}\n")

    # --------------------------------------------------------------------------------------------------------------------
    #                           SAVING GEX DATA INTO MUDATA FILE
    # --------------------------------------------------------------------------------------------------------------------
    print("\n===== SAVING GEX DATA INTO MUDATA FILE =====")

    # If input was h5ad, create a new MuData object with gex as the only modality (for simplicity - not going to happen in scrna pipeline)
    if args.input.suffix == '.h5ad':
        mu = md.MuData({'gex': adata})
    else:
        # If input was h5mu, update the existing mdata with annotated gex
        mu.mod['gex'] = adata

    mu.update()
    print(mu.obs)
    print(mu.var)

    output = str(args.out)
    print("Saving h5mu data to file {}".format(output))
    mu.write(output)
    print(mu)


if __name__ == '__main__':
    main()
