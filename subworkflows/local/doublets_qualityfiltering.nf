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

        if ( !params.skip_doublets ) {
            DOUBLETS (
                ch_convert_concat_filtered
            )
            ch_versions = ch_versions.mix(DOUBLETS.out.versions.first())
            doublets_out = DOUBLETS.out.doublets
            .map { meta, file -> [meta, file] }
            .ifEmpty { [[id: 'dummy'], []] }
        } else {
            doublets_out = [[id: 'dummy'], []]
        }

        //
        // MODULE: Filtered cells of low quality for GEX and CITE modalities in the concatenated h5mu file
        //

        QUALITY_FILTERING (
            ch_h5mu_concat_filtered,
            doublets_out,
            mt_threshold,
        )
        ch_versions = ch_versions.mix(QUALITY_FILTERING.out.versions.first())


    emit:
    ch_versions
    h5mu = QUALITY_FILTERING.out.h5mu



}
