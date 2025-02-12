/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { DOUBLETS } from '../../modules/local/doublets'
include { QUALITY_FILTERING } from '../../modules/local/qualitymetrics'

workflow DOUBLETS_QUALITYFILTERING {

    take:
    ch_convert_concat_filtered
    ch_h5mu_concat_filtered
    mt_threshold

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Compute doublet score for each sample in the concatenated rds file
        //
        DOUBLETS (
            ch_convert_concat_filtered
        )
        ch_versions = ch_versions.mix(DOUBLETS.out.versions.first())
        //
        // MODULE: Filtered cells of low quality in the concatenated h5ad file
        //
        QUALITY_FILTERING (
            ch_h5mu_concat_filtered,
            DOUBLETS.out.doublets,
            mt_threshold
        )
        ch_versions = ch_versions.mix(QUALITY_FILTERING.out.versions.first())
        
    emit:
    ch_versions
    h5mus = QUALITY_FILTERING.out.h5mu

}
