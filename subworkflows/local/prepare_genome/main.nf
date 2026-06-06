/*
 * Prepare reference FASTA and GTF for alignment (gunzip, filter to genome sequences, optional GTF source fix)
 */

include { GUNZIP as GUNZIP_FASTA } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../../../modules/nf-core/gunzip/main'
include { GTF_GENE_FILTER          } from '../../../modules/local/gtf_gene_filter/main'
include { GAWK as GTF_SOURCE_FIX   } from '../../../modules/nf-core/gawk/main'

workflow PREPARE_GENOME {
    take:
    fasta
    gtf
    gtf_source_fix

    main:
    ch_versions     = channel.empty()

    if (fasta) {
        ch_fasta = file(fasta, checkIfExists: true)
        if (fasta.endsWith('.gz')) {
            ch_fasta = GUNZIP_FASTA([[:], ch_fasta]).gunzip
                .map { _meta, fasta_file -> fasta_file }.collect()
        }
        else {
            ch_fasta = channel.value(ch_fasta)
        }
    }

    if (gtf) {
        ch_gtf = file(gtf, checkIfExists: true)
        if (gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF([[:], ch_gtf]).gunzip
                .map { _meta, gtf_file -> gtf_file }.collect()
        }
        else {
            ch_gtf = channel.value(ch_gtf)
        }
    }

    GTF_GENE_FILTER(
        ch_fasta,
        ch_gtf
    )
    ch_gtf = GTF_GENE_FILTER.out.gtf
    ch_versions = ch_versions.mix(GTF_GENE_FILTER.out.versions)

    if (gtf_source_fix) {
        // iGenomes GTF annotations with spaces in the source column (e.g. NCBI GRCh38
        // "Curated Genomic") fail Cell Ranger 10 mkref. Opt-in per genome via
        // gtf_source_has_spaces in the genomes map; see usage docs.
        log.warn(
            "Using an iGenomes GTF with spaces in the source column. nf-core/scrnaseq will rewrite the GTF source field for " +
            "Cell Ranger compatibility, but we recommend current reference annotations for production runs. See " +
            "https://nf-co.re/scrnaseq/dev/docs/usage#reference-genome-options"
        )
        ch_gtf_meta = ch_gtf.map { reference_gtf ->
            [[id: "${reference_gtf.baseName}.source_fixed"], reference_gtf]
        }
        GTF_SOURCE_FIX(ch_gtf_meta, [], false)
        ch_gtf = GTF_SOURCE_FIX.out.output.map { _meta, reference_gtf -> reference_gtf }
        ch_versions = ch_versions.mix(GTF_SOURCE_FIX.out.versions_gawk)
    }

    emit:
    fasta      = ch_fasta
    gtf        = ch_gtf
    versions   = ch_versions
}
