/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { STAR_ALIGN                  } from '../../modules/local/star_align'
include { STAR_GENOMEPARAMS_UPGRADE   } from '../../modules/local/star_genomeparams_upgrade'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { STAR_GENOMEGENERATE }         from '../../modules/nf-core/star/genomegenerate/main'


workflow STARSOLO {
    take:
    genome_fasta
    gtf
    star_index
    star_index_legacy
    protocol
    barcode_whitelist
    ch_fastq
    star_feature
    other_10x_parameters

    main:
    ch_versions = channel.empty()

    assert star_index || (genome_fasta && gtf):
        "Must provide a genome fasta file ('--fasta') and a gtf file ('--gtf') if no index is given!"

    assert gtf: "Must provide a gtf file ('--gtf') for STARSOLO"

    /*
    * Build STAR index if not supplied, or upgrade legacy iGenomes metadata when requested
    */
    if (!star_index) {
        STAR_GENOMEGENERATE(
            genome_fasta.map{ f -> [[id: f.baseName], f]},
            gtf.map{ g -> [[id: g.baseName], g]}
        )
        star_index = STAR_GENOMEGENERATE.out.index.collect()
    }
    else if (star_index_legacy) {
        STAR_GENOMEPARAMS_UPGRADE(star_index)
        ch_versions = ch_versions.mix(STAR_GENOMEPARAMS_UPGRADE.out.versions_gawk)
        star_index = STAR_GENOMEPARAMS_UPGRADE.out.index.collect()
    }

    /*
    * Perform mapping with STAR
    */
    STAR_ALIGN(
        ch_fastq,
        star_index,
        gtf,
        barcode_whitelist,
        protocol,
        star_feature,
        other_10x_parameters
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    raw_counts = STAR_ALIGN.out.raw_counts
        .join(STAR_ALIGN.out.raw_velocyto, remainder: true)
        .map{
            meta, count, velocity ->
                [meta + [input_type: 'raw'], velocity ? [count, velocity] : [count]]
        }

    filtered_counts = STAR_ALIGN.out.filtered_counts
        .join(STAR_ALIGN.out.filtered_velocyto, remainder: true)
        .map{ meta, count, velocity ->
            [meta + [input_type: 'filtered'], velocity ? [count, velocity] : [count]]
        }

    emit:
    ch_versions
    // get rid of meta for star index
    star_result     = STAR_ALIGN.out.tab
    star_counts     = STAR_ALIGN.out.counts
    raw_counts      = raw_counts
    filtered_counts = filtered_counts
    for_multiqc     = STAR_ALIGN.out.log_final.map{ meta, it -> it }
}
