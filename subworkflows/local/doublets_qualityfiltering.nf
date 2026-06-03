/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { DOUBLETS } from '../../modules/local/doublets'
include { QUALITY_FILTERING } from '../../modules/local/qualitymetrics'

workflow DOUBLETS_QUALITYFILTERING {

    take:
    ch_convert_concat_filtered
    ch_h5mu_concat_filtered
    mt_threshold
    min_umi_gex
    max_umi_gex
    min_genes_gex
    max_genes_gex
    min_cells_gex
    min_features_adt
    min_counts_adt


    main:
        ch_versions = channel.empty()

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

        if ( !params.skip_qcfilters ) {
            QUALITY_FILTERING (
                ch_h5mu_concat_filtered,
                doublets_out,
                mt_threshold,
                min_umi_gex,
                max_umi_gex,
                min_genes_gex,
                max_genes_gex,
                min_cells_gex,
                min_features_adt,
                min_counts_adt
            )
            ch_versions = ch_versions.mix(QUALITY_FILTERING.out.versions.first())
            h5mu_out = QUALITY_FILTERING.out.h5mu
        } else {
            h5mu_out = ch_h5mu_concat_filtered
        }

    emit:
    ch_versions
    h5mu = h5mu_out



}
