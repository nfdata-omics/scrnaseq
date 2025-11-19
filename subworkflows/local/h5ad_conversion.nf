/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { CONCAT_H5AD           } from '../../modules/local/concat_h5ad.nf'
include { ANNDATAR_CONVERT      } from '../../modules/local/anndatar_convert'

workflow H5AD_CONVERSION {

    take:
    ch_h5ads
    samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Concat all raw and unfiltered h5ad files
    // Exclude pooled outputs (is_pooled=true) as they're already aggregated by the aligner
    // For aligners that don't set is_pooled, default to including them (is_pooled will be null/false)
    //
    ch_concat_h5ad_input = ch_h5ads
        .filter{ meta, file -> !(meta.is_pooled ?: false) }  // Exclude if is_pooled == true
        .map{ meta, file -> [ [id: 'combined', input_type: meta.input_type], file ]}
        .groupTuple()

    view(ch_concat_h5ad_input)
    CONCAT_H5AD (
        ch_concat_h5ad_input,
        samplesheet
    )

    ch_h5ad_concat = CONCAT_H5AD.out.h5ad
    ch_versions = ch_versions.mix(CONCAT_H5AD.out.versions.first())

    //
    // MODULE: Convert to RDS with AnndataR package
    //
    ANNDATAR_CONVERT (
        ch_h5ads.mix(ch_h5ad_concat)
    )
    ch_versions = ch_versions.mix(ANNDATAR_CONVERT.out.versions.first())

    emit:
    ch_versions
    h5ads = ch_h5ads
}
