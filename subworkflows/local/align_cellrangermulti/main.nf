//
// Include modules
//
include { CELLRANGER_MKGTF                  } from "../../../modules/nf-core/cellranger/mkgtf"
include { CELLRANGER_MKREF                  } from "../../../modules/nf-core/cellranger/mkref"
include { CELLRANGER_MKVDJREF               } from "../../../modules/nf-core/cellranger/mkvdjref"
include { CELLRANGER_MULTI                  } from "../../../modules/nf-core/cellranger/multi"
include { PARSE_CELLRANGERMULTI_SAMPLESHEET } from "../../../modules/local/parse_cellrangermulti_samplesheet"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_MULTI_ALIGN {
    take:
        ch_fasta
        ch_gtf
        ch_fastq
        cellranger_gex_index
        cellranger_vdj_index
        ch_multi_samplesheet

    main:
        //
        // TODO: Include checkers for cellranger multi parameter combinations. For example, when VDJ data is given, require VDJ ref. If FFPE, require frna probe sets, etc.
        //

        // since we merged all data as a meta, now we have a channel per sample, which
        // every item is a meta map for each data-type
        // now we can split it back for passing as input to the module
        ch_fastq
        .flatten()
        .map{ meta ->
            def meta_clone = meta.clone()
            def data_dict  = meta_clone.find{ entry -> entry.key == "${meta_clone.feature_type}" }
            def fastqs = data_dict?.value
            meta_clone.remove( data_dict?.key )
            [ meta_clone, fastqs ]
        }
        .branch {
            meta, fastq ->
                gex: meta.feature_type == "gex"
                    return [ meta, fastq ]
                vdj: meta.feature_type == "vdj"
                    return [ meta, fastq ]
                ab: meta.feature_type == "ab"
                    if ((fastq == file("$projectDir/assets/EMPTY", checkIfExists: true)) || params.fb_reference) { // when empty, should not check for reference
                        return [ meta, fastq ]
                    } else {
                        error ("Antibody reference was not specified. Please provide a reference file for feature barcoding (e.g. antibody measurements).\nPlease refer to https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-feature-ref-csv for more details.")
                    }
                beam: meta.feature_type == "beam"
                    return [ meta, fastq ]
                crispr: meta.feature_type == "crispr"
                    return [ meta, fastq ]
                cmo: meta.feature_type == "cmo"
                    return [ meta, fastq ]
        }
        .set { ch_grouped_fastq }

        // Assign other cellranger reference files
        ch_gex_frna_probeset      = params.gex_frna_probe_set            ? file(params.gex_frna_probe_set)            : []
        ch_gex_target_panel       = params.gex_target_panel              ? file(params.gex_target_panel)              : []
        ch_gex_cmo_set            = params.gex_cmo_set                   ? file(params.gex_cmo_set)                   : []
        ch_gex_barcodes           = params.gex_barcode_sample_assignment ? file(params.gex_barcode_sample_assignment) : []
        ch_fb_reference           = params.fb_reference                  ? file(params.fb_reference)                  : []
        ch_vdj_primer_index       = params.vdj_inner_enrichment_primers  ? file(params.vdj_inner_enrichment_primers)  : []
        ch_beam_antigen_panel_csv = [] // currently not implemented
        ch_beam_control_panel_csv = [] // currently not implemented

        // parse frna and barcode information
        if (ch_multi_samplesheet) {

            //
            // Here, we parse the received cellranger multi barcodes samplesheet.
            // We first use the get the PARSE_CELLRANGERMULTI_SAMPLESHEET module to check it and guarantee structure
            // and also split it to have one fnra/cmo/ocm .csv for each sample.
            //
            // The selection of the GEX fastqs is because samples are always expected to have at least GEX data.
            // Then, using "combined" map, which means, the "additional barcode information" of each sample, we then,
            // parse it to generate the cmo / ocm /frna samplesheets to be used by each sample.
            //
            // Here, to guarantee it and take advantage of the "FIFO"-rule and are sure that the data used in the
            // module is from the same sample from the "normal" samplesheet. We have to use the .concat().groupTuple()
            // pipe instead of .join() because .join() outputs first the arrays that could be joined and afterwards
            // the ones with "remainders", thus, we would not ensure "FIFO" and the same order.
            //
            // To guarantee this, we can define two nf-tests, one having only one sample with CMO and another with two
            // samples using CMOs, even if wrongly/repeated, but just to guarantee FIFO is working.
            //

            PARSE_CELLRANGERMULTI_SAMPLESHEET( ch_multi_samplesheet )

            // CMO
            ch_grouped_fastq.gex
            .map{ pair -> [pair[0].id] }
            .concat( PARSE_CELLRANGERMULTI_SAMPLESHEET.out.cmo.flatten().map { csv -> [ "${csv.baseName}" - "_cmo", csv ] } )
            .groupTuple()
            .map { grp -> if ( grp.size() == 2 ) { grp[1] } else { [] } } // a correct tuple from snippet will have: [ sample, cmo.csv ]
            .set { ch_cmo_barcode_csv }

            // OCM
            ch_grouped_fastq.gex
            .map{ pair -> [pair[0].id] }
            .concat( PARSE_CELLRANGERMULTI_SAMPLESHEET.out.ocm.flatten().map { csv -> [ "${csv.baseName}" - "_ocm", csv ] } )
            .groupTuple()
            .map { grp -> if ( grp.size() == 2 ) { grp[1] } else { [] } } // a correct tuple from snippet will have: [ sample, ocm.csv ]
            .set { ch_ocm_barcode_csv }

            // FRNA
            ch_grouped_fastq.gex
            .map{ pair -> [pair[0].id] }
            .concat( PARSE_CELLRANGERMULTI_SAMPLESHEET.out.frna.flatten().map { csv -> [ "${csv.baseName}" - "_frna", csv ] } )
            .groupTuple()
            .map { grp -> if ( grp.size() == 2 ) { grp[1] } else { [] } } // a correct tuple from snippet will have: [ sample, frna.csv ]
            .set { ch_frna_sample_csv }

            ch_grouped_fastq.gex.view()
            PARSE_CELLRANGERMULTI_SAMPLESHEET.out.frna.flatten().view()
            ch_frna_sample_csv.view()

        } else {
            ch_cmo_barcode_csv = []
            ch_ocm_barcode_csv = []
            ch_frna_sample_csv = []
        }

        //
        // Prepare GTF
        //
        if ( !cellranger_gex_index || (!cellranger_vdj_index && !params.skip_cellrangermulti_vdjref) ) {

            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_MKGTF ( ch_gtf )

        }

        //
        // Prepare gex reference (Normal Ref)
        //
        if ( !cellranger_gex_index ) {

            // Validate that gex_reference_version is provided when required
            if ( params.gex_frna_probe_set && !params.gex_reference_version ) {
                error "Parameter 'gex_reference_version' is required when 'gex_frna_probe_set' is provided and 'cellranger_index' is not provided. The reference genome version must match the probeset reference."
            }

            // Validate that gex_reference_version matches the probeset reference genome
            if ( params.gex_frna_probe_set && params.gex_reference_version ) {
                def probeset_file = file(params.gex_frna_probe_set)
                def probeset_reference = null
                def done = false
                probeset_file.eachLine { line ->
                    if (done)
                        return
                    if (line.startsWith("#reference_genome=")) {
                        def ref_split = line.split("=")
                        if (ref_split.size() > 1) {
                            probeset_reference = ref_split[1].trim()
                        }
                    }
                }
                if ( probeset_reference != params.gex_reference_version ) {
                    error "Parameter 'gex_reference_version' (${params.gex_reference_version}) does not match the probeset reference genome (${probeset_reference}). Please ensure the reference genome version matches the probeset file."
                }
            }

            // Make reference genome
            def reference_name = params.gex_reference_version ?: "gex_reference_version"
            CELLRANGER_MKREF(
                ch_fasta,
                CELLRANGER_MKGTF.out.gtf,
                reference_name
            )
            ch_cellranger_gex_index = CELLRANGER_MKREF.out.reference.ifEmpty { [] }

        } else {
            ch_cellranger_gex_index = cellranger_gex_index
        }

        //
        // Prepare vdj reference (Special)
        //
        if ( !cellranger_vdj_index ) {

            if ( !params.skip_cellrangermulti_vdjref  ) { // if user uses cellranger multi but does not have VDJ data
                // Make reference genome
                CELLRANGER_MKVDJREF(
                    ch_fasta,
                    CELLRANGER_MKGTF.out.gtf,
                    [], // currently ignoring the 'seqs' option
                    "vdj_reference"
                )
                ch_cellranger_vdj_index = CELLRANGER_MKVDJREF.out.reference.ifEmpty { [] }
            } else {
                ch_cellranger_vdj_index = []
            }

        } else {
            ch_cellranger_vdj_index = cellranger_vdj_index
        }

        //
        // MODULE: cellranger multi
        //
        CELLRANGER_MULTI(
            ch_grouped_fastq.gex.map{ pair -> pair[0] },
            ch_grouped_fastq.gex.map { meta, fastqs -> [meta, fastqs, meta.options] },
            ch_grouped_fastq.vdj,
            ch_grouped_fastq.ab,
            ch_grouped_fastq.beam,
            ch_grouped_fastq.cmo,
            ch_grouped_fastq.crispr,
            ch_cellranger_gex_index,
            ch_gex_frna_probeset,
            ch_gex_target_panel,
            ch_cellranger_vdj_index,
            ch_vdj_primer_index,
            ch_fb_reference,
            ch_beam_antigen_panel_csv,
            ch_beam_control_panel_csv,
            ch_gex_cmo_set,
            ch_cmo_barcode_csv,
            [],
            ch_frna_sample_csv,
            ch_ocm_barcode_csv,
            params.skip_cellranger_renaming
        )

        //
        // Cell Ranger multi produces aggregate matrices for the complete run and sample-specific matrices under
        // per_sample_outs. Only the latter belong in the downstream sample channel: aggregate matrices do not
        // represent biological samples and would duplicate cells already present in the per-sample outputs.
        //
        ch_matrices_filtered = parse_per_sample_output_channels(
            CELLRANGER_MULTI.out.outs,
            "sample_filtered_feature_bc_matrix"
        )
        ch_matrices_raw = parse_per_sample_output_channels(
            CELLRANGER_MULTI.out.outs,
            "sample_raw_feature_bc_matrix"
        )

        // Preserve the sample and receptor type associated with each filtered V(D)J annotation.
        ch_vdj_files = CELLRANGER_MULTI.out.outs
            .flatMap { meta, outs ->
                outs.findAll { path ->
                    path.name == "filtered_contig_annotations.csv" &&
                    path.toString().contains("/per_sample_outs/") &&
                    (path.parent.name in ["vdj_b", "vdj_t"])
                }.collect { path ->
                    def meta_clone = per_sample_meta(meta, path)
                    meta_clone.feature_type = "vdj"
                    meta_clone.vdj_type = path.parent.name
                    [meta_clone, path]
                }
            }
            .collect()
            .map { entries ->
                def pairs = entries.collate(2)
                [
                    pairs.collect { meta, _path -> meta },
                    pairs.collect { _meta, path -> path }
                ]
            }

    emit:
        cellrangermulti_out          = CELLRANGER_MULTI.out.outs
        cellrangermulti_mtx_raw      = ch_matrices_raw
        cellrangermulti_mtx_filtered = ch_matrices_filtered
        vdj                          = ch_vdj_files
}

def parse_per_sample_output_channels(in_ch, matrix_name) {
    return in_ch
        .flatMap { meta, matrix_files ->
            matrix_files.findAll { path ->
                path.toString().contains("/per_sample_outs/") &&
                (path.name == matrix_name || path.name == "${matrix_name}.h5")
            }.collect { path ->
                def meta_clone = per_sample_meta(meta, path)
                meta_clone.feature_type = "gex"
                meta_clone.input_type = matrix_name.contains("_raw_") ? "raw" : "filtered"
                [meta_clone, path]
            }
        }
        .groupTuple(by: 0)
}

def per_sample_meta(meta, path) {
    def meta_clone = meta.clone()
    meta_clone.id = path.toString().split("/per_sample_outs/")[1].split("/")[0]
    return meta_clone
}
