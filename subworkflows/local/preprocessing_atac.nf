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
    cell_annotation_meta_ch

    

    main:
        ch_versions = Channel.empty()

        
        QUALITY_FILTERING_ATAC (
            fragments,
            fragments_index,
            nucleosome_threshold,
            tss_threshold,
            blacklist_path

        )
        ch_versions = ch_versions.mix(QUALITY_FILTERING_ATAC.out.versions)
        
    
        DIMENSIONALITY_REDUCTION_ATAC (
            QUALITY_FILTERING_ATAC.out.h5ad
        )
        ch_versions = ch_versions.mix(DIMENSIONALITY_REDUCTION_ATAC.out.versions)

        
        PEAK_CALLING(
            DIMENSIONALITY_REDUCTION_ATAC.out.h5ad,
            cell_annotation_meta_ch
        )
        ch_versions = ch_versions.mix(PEAK_CALLING.out.versions)
        

    emit:
    ch_versions
    h5ad = PEAK_CALLING.out.h5ad_peak
    


}
