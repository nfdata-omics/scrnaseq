/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { QUALITY_FILTERING_ATAC } from '../../modules/local/qualitymetrics_atac'
include { DOUBLETS_ATAC } from '../../modules/local/doublets_atac'

workflow ATAC_PREPROCESSING {

    take:
    fragments
    fragments_index
    nucleosome_threshold
    tss_threshold
    

    main:
        ch_versions = Channel.empty()

        
        QUALITY_FILTERING_ATAC (
            fragments,
            fragments_index,
            nucleosome_threshold,
            tss_threshold
        )
        ch_versions = ch_versions.mix(QUALITY_FILTERING_ATAC.out.versions.first())
        
        DOUBLETS_ATAC (
            QUALITY_FILTERING_ATAC.out.h5ad
        )
        ch_versions = ch_versions.mix(DOUBLETS_ATAC.out.versions.first())

        '''
        QUALITY_FILTERING_ATAC (
            fragments,
            fragments_index,
            nucleosome_threshold,
            tss_threshold
        )
        ch_versions = ch_versions.mix(QUALITY_FILTERING_ATAC.out.versions.first())
        '''

    emit:
    ch_versions
    h5ad = DOUBLETS_ATAC.out.h5ad


}
