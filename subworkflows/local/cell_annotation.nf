include { CELLANNOTATION_CELLTYPIST       } from '../../modules/local/cellannotation/main.nf'
include { CELL_ANNOTATION_LLM              } from '../../modules/local/cellannotation_llm/main.nf'

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

    CELLANNOTATION_LLM(
        h5mu,
        resolutions,
        llm_species,
        llm_tissue
    )

    ch_versions = ch_versions.mix(CELLANNOTATION_LLM.out.versions)
    ch_llm_h5mu = CELLANNOTATION_LLM.out.h5mu
    ch_llm_umap = CELLANNOTATION_LLM.out.graph_umap
    ch_llm_annotated_clusters = CELLANNOTATION_LLM.out.annotated_clusters
    ch_llm_parameters_txt = CELLANNOTATION_LLM.out.parameters_txt

    emit:

    ch_versions

    ch_mu5ad
    cell_annotation_meta_ch
    ch_versions
    ch_llm_h5mu
    ch_llm_umap
    ch_llm_annotated_clusters
    ch_llm_parameters_txt

}
