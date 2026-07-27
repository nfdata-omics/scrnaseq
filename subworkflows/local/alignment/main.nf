include { CELLRANGER_ALIGN                                  } from "../align_cellranger"
include { CELLRANGER_MULTI_ALIGN                            } from "../align_cellrangermulti"
include { CELLRANGERARC_ALIGN                               } from "../align_cellrangerarc"
include { cellrangerarcStructure                            } from '../utils_nfcore_scrnaseq_pipeline'

workflow ALIGNMENT {
    take:
        ch_samplesheet
        ch_genome_fasta
        ch_genome_gtf
        ch_cellranger_index
        protocol
        ch_motifs
        ch_cellrangerarc_config
        cellranger_vdj_index
        ch_multi_samplesheet

    main:

        ch_multiqc_files = channel.empty()
        ch_vdj_file      = channel.empty()
        ch_transformed_fragments_channel = channel.empty()
        ch_transformed_fragments_index_channel = channel.empty()

        empty_file = file("$projectDir/assets/EMPTY", checkIfExists: true)

        // filter the samplesheet to obtained files that have been already pre-processed
        ch_samplesheet
            .filter { meta, _files -> meta.input_type in ['raw', 'filtered'] }
            .set { ch_mtx_matrices }

        // extract the fastq files from the samplesheet to be used in the alignment subworkflows
        ch_fastq = ch_samplesheet.filter { meta, _files -> meta.input_type == 'fastq' }

        if (params.aligner == "cellranger") {

            CELLRANGER_ALIGN(
                ch_genome_fasta,
                ch_genome_gtf,
                ch_cellranger_index,
                ch_fastq,
                protocol
            )
            ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGER_ALIGN.out.cellranger_matrices_raw, CELLRANGER_ALIGN.out.cellranger_matrices_filtered )
            ch_multiqc_files = ch_multiqc_files.mix(CELLRANGER_ALIGN.out.cellranger_out.map {
                _meta, outs -> outs.findAll{ summary -> summary.name == "web_summary.html"}
            })
        }

        if (params.aligner == "cellrangerarc") {

            ch_cellrangerarc_fastq = ch_fastq
                .flatMap { meta, fastqs ->
                    def library_size = meta.feature_type == 'atac' ? 3 : 2
                    if (fastqs.size() % library_size != 0) {
                        error("Please check input samplesheet -> Unexpected number of FASTQ files for ${meta.id} (${meta.feature_type}).")
                    }
                    fastqs.collate(library_size).collect { library_fastqs ->
                        [meta.id, meta, library_fastqs]
                    }
                }
                .groupTuple()
                .map { grouped_fastqs ->
                    cellrangerarcStructure(grouped_fastqs)
                }

            CELLRANGERARC_ALIGN(
                ch_genome_fasta,
                ch_genome_gtf,
                ch_motifs,
                ch_cellranger_index,
                ch_cellrangerarc_fastq,
                ch_cellrangerarc_config
            )
            ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGERARC_ALIGN.out.cellrangerarc_mtx_raw, CELLRANGERARC_ALIGN.out.cellrangerarc_mtx_filtered )

            // Collect the fragments files and their index
            ch_fragments =
                CELLRANGERARC_ALIGN.out.cellrangerarc_out.map { meta, outs ->
                    def desired_files = outs.findAll { file -> file.name == "atac_fragments.tsv.gz" }
                    if (desired_files.size() > 0) {
                        [meta, desired_files]
                    }
                    else {
                    }
            }
            ch_fragments_collect =  ch_fragments.collect()

            ch_transformed_fragments_channel = ch_fragments_collect.map { list ->
                def meta = []
                def files = []

                list.collate(2).each { pair ->
                    meta << pair[0]
                    files << pair[1]
                }
                return [meta, files.flatten()]
            }

            ch_fragments_index = CELLRANGERARC_ALIGN.out.cellrangerarc_out.map { meta, outs ->
                def desired_files = outs.findAll { file -> file.name == "atac_fragments.tsv.gz.tbi" }

                if (desired_files.size() > 0) {
                    [meta, desired_files]
                }
                else {
                }
            }
            ch_vdj_fragments_index_collect =  ch_fragments_index.collect()

            ch_transformed_fragments_index_channel = ch_vdj_fragments_index_collect.map { list ->
                def meta = []
                def files = []

                list.collate(2).each { pair ->
                    meta << pair[0]
                    files << pair[1]
                }
                return [meta, files.flatten()]
                }
        }

        // Run cellrangermulti pipeline
        if (params.aligner == 'cellrangermulti') {

            // parse the input data to generate a collected channel per sample, which will have
            // the metadata and data for each data-type of every sample.
            // then, inside the subworkflow, it can be parsed to manage inputs to the module
            ch_fastq
            .map { meta, fastqs ->
                def parsed_meta = meta.clone() + [ "${meta.feature_type.toString()}": fastqs ]
                parsed_meta.options = [:]

                // add an universal key to differentiate from empty channels so that the "&& meta_gex?.options" lines in the module main.nf can work properly
                parsed_meta.options['data-available'] = true

                // add cellranger options that are currently handled by pipeline, coming from samplesheet
                // the module parses them from the 'gex' options
                if (meta.feature_type.toString() == 'gex') {
                    parsed_meta.options['create-bam'] = params.save_align_intermeds  // force bam creation -- param required by cellranger multi
                    if (meta.expected_cells) { parsed_meta.options['expected-cells'] = meta.expected_cells }
                    parsed_meta.options['chemistry'] = protocol
                }

                [ parsed_meta.id , parsed_meta ]
            }
            .groupTuple( by: 0 )
            .map{ sample_id, map_collection ->
                // Now we must check if every data possibility taken into account in the .branch() operation
                // performed inside the CELLRANGER_MULTI_ALIGN subworkflow are initialized, even with empty files
                // This to ensure that the sizes of each data channel is the same, and the the order and the data types
                // are used together with its rightful pairs
                //
                // data.types: gex, vdj, ab, beam, crispr, cmo

                // clone ArrayBag (received from .groupTuple()) to avoid mutating the input
                def map_collection_clone = []
                map_collection_clone.addAll(map_collection)

                // generate the expected EMPTY tuple when a data type is not used
                // needs to have a collected map like that, so every sample from the samplesheet is analysed one at a time,
                // allowing to have multiple samples in the sheet, having all the data-type tuples initialized,
                // either empty or populated. It will be branched inside the subworkflow.
                if (!map_collection_clone.any{ m -> m.feature_type == 'gex' })    { map_collection_clone.add( [id: sample_id, feature_type: 'gex'   , gex:    empty_file, options:[:] ] ) }
                if (!map_collection_clone.any{ m -> m.feature_type == 'vdj' })    { map_collection_clone.add( [id: sample_id, feature_type: 'vdj'   , vdj:    empty_file, options:[:] ] ) }
                if (!map_collection_clone.any{ m -> m.feature_type == 'ab' })     { map_collection_clone.add( [id: sample_id, feature_type: 'ab'    , ab:     empty_file, options:[:] ] ) }
                if (!map_collection_clone.any{ m -> m.feature_type == 'beam' })   { map_collection_clone.add( [id: sample_id, feature_type: 'beam'  , beam:   empty_file, options:[:] ] ) } // currently not implemented, the input samplesheet checking will not allow it.
                if (!map_collection_clone.any{ m -> m.feature_type == 'crispr' }) { map_collection_clone.add( [id: sample_id, feature_type: 'crispr', crispr: empty_file, options:[:] ] ) }
                if (!map_collection_clone.any{ m -> m.feature_type == 'cmo' })    { map_collection_clone.add( [id: sample_id, feature_type: 'cmo'   , cmo:    empty_file, options:[:] ] ) }

                // return final map
                map_collection_clone
            }
            .set{ ch_cellrangermulti_collected_channel }

            // Run cellranger multi
            CELLRANGER_MULTI_ALIGN(
                ch_genome_fasta,
                ch_genome_gtf,
                ch_cellrangermulti_collected_channel,
                //ch_transformed_fragments_index_channel,
                ch_cellranger_index,
                cellranger_vdj_index,
                ch_multi_samplesheet
            )
            ch_multiqc_files = ch_multiqc_files.mix( CELLRANGER_MULTI_ALIGN.out.cellrangermulti_out.map{
                _meta, outs -> outs.findAll{ it -> it.name == "web_summary.html" }
            })
            ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGER_MULTI_ALIGN.out.cellrangermulti_mtx_raw, CELLRANGER_MULTI_ALIGN.out.cellrangermulti_mtx_filtered )
            ch_vdj_file     = CELLRANGER_MULTI_ALIGN.out.vdj

        }

    emit:
        mtx_matrices = ch_mtx_matrices
        multiqc_files = ch_multiqc_files
        vdj_file = ch_vdj_file
        fragments_file = ch_transformed_fragments_channel
        fragments_index = ch_transformed_fragments_index_channel

}
