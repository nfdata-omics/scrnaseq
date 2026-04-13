include { CELLANNOTATION_CELLTYPIST       } from '../../modules/local/cellannotation/main.nf'
include { CELL_ANNOTATION_LLM             } from '../../modules/local/cellannotation_llm/main.nf'
include { CUSTOM_GENES                    } from '../../modules/local/custom_genes'

workflow CELL_ANNOTATION {

    take:

    h5mu
    input_model
    resolutions
    llm_species
    llm_tissue

    main:

    ch_versions = Channel.empty()

    if ( params.input_model ) {
        CELLANNOTATION_CELLTYPIST (
            h5mu,
            input_model,
            params.annotation_method,
            params.annotation_threshold
        )
        ch_versions = ch_versions.mix(CELLANNOTATION_CELLTYPIST.out.versions)
        ch_mu5ad = CELLANNOTATION_CELLTYPIST.out.h5mu
        cell_annotation_meta_ch = CELLANNOTATION_CELLTYPIST.out.metadata
    } else {
        ch_mu5ad = h5mu
        cell_annotation_meta_ch = channel.empty()
    }

    if ( params.llm_species && params.llm_tissue) {

        CELL_ANNOTATION_LLM(
            h5mu,
            resolutions,
            llm_species,
            llm_tissue
        )

        ch_versions = ch_versions.mix(CELL_ANNOTATION_LLM.out.versions)
        ch_llm_h5mu = CELL_ANNOTATION_LLM.out.h5mu
        ch_llm_umap = CELL_ANNOTATION_LLM.out.graph_umap
        ch_llm_annotated_clusters = CELL_ANNOTATION_LLM.out.annotated_clusters
        ch_llm_parameters_txt = CELL_ANNOTATION_LLM.out.parameters_txt

    } else {
        ch_llm_h5mu = channel.empty()
        ch_llm_umap = channel.empty()
        ch_llm_annotated_clusters = channel.empty()
        ch_llm_parameters_txt = channel.empty()
    }

    if ( params.custom_geneset ) {
        ch_custom_geneset = Channel.fromList(params.custom_geneset.split(',').flatten())

        if ( params.resolution ) {
            resolution_ch = Channel.fromList(params.resolution.toString().split(',').flatten())

            resolution_ch
                .combine( ch_custom_geneset )
                .map{ res, genes -> [["res": res, "genes": genes], res, genes] }
                .set { ch_res_geneset }
        } else {
            // if no resolution is provided, use 100 as fake resolution
            fake_res = 100
            ch_res_geneset = ch_custom_geneset.map { genes ->
                [["res": fake_res, "genes": genes], fake_res, genes]
            }
        }
        CUSTOM_GENES (
            h5mu,
            ch_res_geneset
        )
        ch_versions = ch_versions.mix(CUSTOM_GENES.out.versions)

        ch_featplot = CUSTOM_GENES.out.feat_plot
        ch_dotplot = CUSTOM_GENES.out.dotplot
        ch_heatmap = CUSTOM_GENES.out.heatmap
        ch_violin = CUSTOM_GENES.out.violin
    }

    emit:

    ch_versions
    ch_mu5ad
    cell_annotation_meta_ch
    ch_llm_h5mu
    ch_llm_umap
    ch_llm_annotated_clusters
    ch_llm_parameters_txt
    ch_featplot
    ch_dotplot
    ch_heatmap
    ch_violin

}
