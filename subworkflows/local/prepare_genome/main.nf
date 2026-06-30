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
    ch_versions = channel.empty()
    ch_fasta    = []
    ch_gtf      = []

    if (fasta) {
        def fasta_file = file(fasta, checkIfExists: true)
        if (fasta.endsWith('.gz')) {
            ch_fasta = GUNZIP_FASTA([[:], fasta_file]).gunzip
                .map { _meta, fasta_path -> fasta_path }.collect()
        }
        else {
            ch_fasta = channel.value(fasta_file)
        }
    }

    if (gtf) {
        def gtf_file = file(gtf, checkIfExists: true)
        if (gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF([[:], gtf_file]).gunzip
                .map { _meta, gtf_path -> gtf_path }.collect()
        }
        else {
            ch_gtf = channel.value(gtf_file)
        }

        if (fasta) {
            GTF_GENE_FILTER(
                ch_fasta,
                ch_gtf
            )
            ch_gtf = GTF_GENE_FILTER.out.gtf
            ch_versions = ch_versions.mix(GTF_GENE_FILTER.out.versions)
        }

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
        }
    }

    emit:
    fasta      = ch_fasta
    gtf        = ch_gtf
    versions   = ch_versions
}
