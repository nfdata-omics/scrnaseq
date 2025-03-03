/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { NORMALIZATION } from '../../modules/local/normalization'
include { HIGHLY_VARIABLE_GENES } from '../../modules/local/highly_variable_genes'

workflow NORMALIZATION_AND_HVG {

    take:
    h5mus
    ch_h5ad_concat_raw

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Normalize count matrices contained in the concatenated h5ad file
        //
        NORMALIZATION (
            h5mus,
            ch_h5ad_concat_raw
        )
        ch_versions = ch_versions.mix(NORMALIZATION.out.versions.first())

        //
        // MODULE: Highly variable genes detection, added with gene annotation
        //
        HIGHLY_VARIABLE_GENES (
            NORMALIZATION.out.h5mu
        )
        ch_versions = ch_versions.mix(HIGHLY_VARIABLE_GENES.out.versions.first())

    emit:
    ch_versions
    h5mus = HIGHLY_VARIABLE_GENES.out.h5mu

}
