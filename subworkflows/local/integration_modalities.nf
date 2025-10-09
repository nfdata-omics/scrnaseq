/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { INTEGRATION } from '../../modules/local/integration'
include { MOFA_INTEGRATION } from '../../modules/local/MOFA_integration'

workflow INTEGRATION_MODALITIES {

    take:
    h5mu
    h5ad


    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Integrate GEX and ADT modalities indipendently with Harmony
        //
        INTEGRATION (
            h5mu
        )
        ch_versions = ch_versions.mix(INTEGRATION.out.versions.first())

        //
        // MODULE: Compute WNN graph for each modality
        //
        if (!params.skip_integration){
            MOFA_INTEGRATION (
                INTEGRATION.out.h5mu,
                h5ad
            )
            ch_versions = ch_versions.mix(MOFA_INTEGRATION.out.versions.first())

            emit:
            ch_versions
            h5mu = MOFA_INTEGRATION.out.h5mu

         } else {
            emit:
            ch_versions
            h5mu = INTEGRATION.out.h5mu
         }

        emit:
        ch_versions
        h5mu

}
