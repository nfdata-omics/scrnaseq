/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { QUALITY_FILTERING_ATAC } from '../../modules/local/qualitymetrics_atac'
include { DIMENSIONALITY_REDUCTION_ATAC } from '../../modules/local/dimensionality_reduction_atac'
include { PEAK_CALLING } from '../../modules/local/peak_calling'

workflow ATAC_PREPROCESSING {

    take:
    fragments
    fragments_index
    nucleosome_threshold
    tss_threshold
    blacklist_path
    

    main:
        ch_versions = Channel.empty()

        
        QUALITY_FILTERING_ATAC (
            fragments,
            fragments_index,
            nucleosome_threshold,
            tss_threshold,
            blacklist_path

        )
        ch_versions = ch_versions.mix(QUALITY_FILTERING_ATAC.out.versions.first())
        
    
        DIMENSIONALITY_REDUCTION_ATAC (
            QUALITY_FILTERING_ATAC.out.h5ad
        )
        ch_versions = ch_versions.mix(DIMENSIONALITY_REDUCTION_ATAC.out.versions.first())

        '''
        PEAK_CALLING(
            DIMENSIONALITY_REDUCTION_ATAC.out.h5ad
        )
        ch_versions = ch_versions.mix(PEAK_CALLING.out.versions.first())
        '''

    emit:
    ch_versions
    h5ad = DIMENSIONALITY_REDUCTION_ATAC.out.h5ad


}
