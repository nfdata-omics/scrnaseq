/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { INTEGRATION } from '../../modules/local/integration'
include { WNN_INTEGRATION } from '../../modules/local/WNN_integration'

workflow INTEGRATION_MODALITIES {

    take:
    h5mu

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
        WNN_INTEGRATION (
            INTEGRATION.out.h5mu
        )
        ch_versions = ch_versions.mix(WNN_INTEGRATION.out.versions.first())

    emit:
    ch_versions
    h5mu = WNN_INTEGRATION.out.h5mu

}
