/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { getProtocol                                       } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { gtfSourceFixNeeded                                } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { PREPARE_GENOME                                    } from '../subworkflows/local/prepare_genome'
include { FASTQC_CHECK                                      } from '../subworkflows/local/fastqc'
include { CELLRANGER_ALIGN                                  } from "../subworkflows/local/align_cellranger"
include { CELLRANGER_MULTI_ALIGN                            } from "../subworkflows/local/align_cellrangermulti"
include { CELLRANGERARC_ALIGN                               } from "../subworkflows/local/align_cellrangerarc"
include { MTX_TO_H5AD                                       } from '../modules/local/mtx_to_h5ad'
include { H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA } from '../subworkflows/nf-core/h5ad_removebackground_barcodes_cellbender_anndata'
include { H5AD_CONVERSION                                   } from '../subworkflows/local/h5ad_conversion'
include { ATAC_PREPROCESSING                                } from '../subworkflows/local/preprocessing_atac'
include { CONCATENATE_VDJ                                   } from '../modules/local/concatenate_vdj'
include { CONVERT_MUDATA                                    } from '../modules/local/convert_mudata'
include { DOUBLETS_QUALITYFILTERING                         } from '../subworkflows/local/doublets_qualityfiltering'
include { NORMALIZATION_AND_HVG                             } from '../subworkflows/local/normalization_and_hvg'
include { CELL_ANNOTATION                                   } from '../modules/local/cellannotation'
include { INTEGRATION_MODALITIES                            } from '../subworkflows/local/integration_modalities'
include { CLUSTERING                                        } from '../modules/local/clustering'
include { CLUSTREE                                          } from '../modules/local/clustree'
include { ENRICH_MARKERS                                    } from '../modules/local/enrich_markers'
include { CUSTOM_GENES                                      } from '../modules/local/custom_genes'
include { DIFFERENTIAL_ABUNDANCE                            } from '../modules/local/differential_abundance'
include { PSEUDOBULK_ANALYSIS                               } from '../subworkflows/local/pseudobulk_analysis'
include { CELL_INTERACTION                                  } from '../modules/local/cell_interaction'

workflow SCRNASEQ {

    take:
    ch_fastq                    // channel: [ meta, fastq ] from samplesheet
    counts
    h5ad_matrix
    fasta                       // val: path-like string (or null)
    gtf                         // val: path-like string (or null)
    gff                         // val: path-like string (or null)
    cellranger_index            // val: path-like string (or null)
    txp2gene                    // val: path-like string (or null)
    motifs                      // val: path-like string (or null)
    cellranger_vdj_index        // val: path-like string (or null)
    multiqc_config              // val: path-like string (or null)
    multiqc_logo                // val: path-like string (or null)
    multiqc_methods_description // val: path-like string (or null)
    outdir                      // val: string

    main:
    ch_multiqc_files = channel.empty()
    ch_versions      = channel.empty()
    ch_mtx_matrices  = channel.empty()

    protocol_config = getProtocol(workflow, log, params.aligner, params.protocol)
    if (protocol_config['protocol'] == 'auto' && params.aligner !in ["cellranger", "cellrangerarc", "cellrangermulti"]) {
        error "Only cellranger supports `protocol = 'auto'`. Please specify the protocol manually!"
    }

    // general input and params
    ch_motifs               = motifs           ? file(motifs, checkIfExists: true)           : []
    ch_txp2gene             = txp2gene         ? file(txp2gene, checkIfExists: true)         : []

    // Warn if both GTF and GFF files are provided
    if (gtf && gff) {
        log.warn("Both GTF and GFF files are provided. GTF file will be used.")
    }

    // samplesheet - this is passed to the MTX conversion functions to add metadata to the
    // AnnData objects.
    ch_input = params.input                ? file(params.input, checkIfExists: true)    : []
    ch_counts = params.counts              ? file(params.counts, checkIfExists: true)    : []
    ch_h5ad_matrix = params.h5ad_matrix    ? file(params.h5ad_matrix, checkIfExists: true): []

    //cellranger params
    ch_cellranger_index = cellranger_index ? file(cellranger_index, checkIfExists: true) : []

    //cellrangermulti params
    cellranger_vdj_index = cellranger_vdj_index             ? file(cellranger_vdj_index, checkIfExists: true)             : []
    ch_multi_samplesheet = params.cellranger_multi_barcodes ? file(params.cellranger_multi_barcodes, checkIfExists: true) : []
    empty_file           = file("$projectDir/assets/EMPTY", checkIfExists: true)

    // cellrangerarc params
    ch_cellrangerarc_config = params.cellrangerarc_config ? file(params.cellrangerarc_config)          : []

    // Differential analysis params
    ch_diff_abundance_comparisons = params.diff_abundance_comparisons ?
        channel.fromList(params.diff_abundance_comparisons.split(',').flatten()) : channel.empty()

    // Pseudobulk params
    ch_pseudobulk_group = params.pseudobulk_group ?
        channel.value(params.pseudobulk_group) : channel.empty()
    ch_pseudobulk_comparisons = params.pseudobulk_comparisons ?
        channel.fromList(params.pseudobulk_comparisons.split(',').flatten()) : channel.empty()
    ch_pseudobulk_formula = params.pseudobulk_formula ?
        channel.value(params.pseudobulk_formula) : channel.empty()
    ch_pseudobulk_fdr = params.pseudobulk_fdr ?
        channel.value(params.pseudobulk_fdr) : channel.empty()

    // Cell interaction params
    ch_liana_method = params.liana_method ?
        channel.value(params.liana_method) : channel.empty()
    ch_liana_resource = params.liana_resource ?
        channel.value(params.liana_resource) : channel.empty()

    // Run FastQC
    if (!params.skip_fastqc) {
        FASTQC_CHECK ( ch_fastq )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CHECK.out.fastqc_multiqc.flatten())
    }

    //
    // Prepare reference FASTA and GTF (gunzip, filter, optional Cell Ranger GTF source fix)
    //
    PREPARE_GENOME(
        fasta,
        gtf,
        gff,
        gtfSourceFixNeeded(params.aligner, params.genome, params.genomes, gtf)
    )
    ch_genome_fasta = PREPARE_GENOME.out.fasta
    ch_genome_gtf   = PREPARE_GENOME.out.gtf

    // Run cellranger pipeline
    if (params.aligner == "cellranger") {
        CELLRANGER_ALIGN(
            ch_genome_fasta,
            ch_genome_gtf,
            ch_cellranger_index,
            ch_fastq,
            protocol_config['protocol']
        )
        ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGER_ALIGN.out.cellranger_matrices_raw, CELLRANGER_ALIGN.out.cellranger_matrices_filtered )
        ch_multiqc_files = ch_multiqc_files.mix(CELLRANGER_ALIGN.out.cellranger_out.map {
            _meta, outs -> outs.findAll{ summary -> summary.name == "web_summary.html"}
        })
    }

    // Run cellrangerarc pipeline
    if (params.aligner == "cellrangerarc") {
        CELLRANGERARC_ALIGN(
            ch_genome_fasta,
            ch_genome_gtf,
            ch_motifs,
            ch_cellranger_index,
            ch_fastq,
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


        ch_fragments_index =
            CELLRANGERARC_ALIGN.out.cellrangerarc_out.map { meta, outs ->
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
                parsed_meta.options['chemistry'] = protocol_config['protocol']
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

    }

    ch_count_matrix = channel.empty()
    if ( params.counts ) {
        ch_count_matrix = channel
        .fromPath(params.counts, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id         : row.sample,
                input_type : row.input_type
            ]
            def matrix_file = file(row.h5)
            tuple(meta, matrix_file)
        }
    } else {
        ch_count_matrix = ch_mtx_matrices
    }


    //
    // MODULE: Convert mtx matrices to h5ad
    //
    MTX_TO_H5AD (
        ch_count_matrix,
        ch_txp2gene,
        [],
        params.aligner
    )
    ch_versions = ch_versions.mix(MTX_TO_H5AD.out.versions.first())
    ch_h5ads = MTX_TO_H5AD.out.h5ad


    //
    // SUBWORKFLOW: Run cellbender remove background subworkflow
    //
    if ( !params.skip_cellbender && !(params.aligner in ['cellrangerarc']) ) {
        // module should only run on the raw matrices thus, filter-out the filtered result of the aligners that can produce it
        H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA (
            ch_h5ads
                .filter { meta, _mtx_files -> meta.input_type == 'raw' }
                .map { meta, mtx_files -> [ meta + [input_type: 'cellbender_filter'], mtx_files ]} // to avoid name collision
        )
        ch_h5ads = ch_h5ads.mix(
            H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA.out.h5ad
        )
    }

    //
    // SUBWORKFLOW: Concat samples and convert h5ad to other formats
    //
    H5AD_CONVERSION (
        ch_h5ads,
        ch_input ?: ch_counts
    )
    ch_versions = ch_versions.mix(H5AD_CONVERSION.out.ch_versions)

    //
    // MODULE: Concat vdj samples and save as h5ad format
    //

    if (params.aligner == "cellrangermulti") {
        CONCATENATE_VDJ (
            CELLRANGER_MULTI_ALIGN.out.vdj
        )
        ch_versions = ch_versions.mix(CONCATENATE_VDJ.out.versions)


    //
    // SUBWORKFLOW: Concat GEX, VDJ and CITE data and save as MuData object
    //
        ch_vdj = CONCATENATE_VDJ.out.h5ad
            .map { meta, file -> [meta, file] }
            .ifEmpty { [[id: 'dummy'], []] }
    } else {
        ch_vdj = [[id: 'dummy'], []]
    }

    //TODO: modify this part beacuse only one input is present
    if (params.demultiplexing_doublets) {
    ch_metadata_demuxafy = channel.fromPath(params.demultiplexing_doublets, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [ id: row.sample ]
            def metadata_file = file(row.path)
            tuple(meta, metadata_file)
        }
    } else {
        ch_metadata_demuxafy = channel.value([ [id: 'dummy'], [] ])
    }

    ch_metadata = params.metadata ? channel.value(params.metadata) : channel.value(file('dummy_metadata.csv'))


    if (params.aligner == "cellrangermulti" || params.aligner == "cellrangerarc" || params.aligner == "cellranger" ) {
        def ch_h5ad_selected = params.counts ?
            H5AD_CONVERSION.out.h5ad_cellbender :
            (
                params.h5ad_matrix ?
                    channel
                        .fromPath(params.h5ad_matrix,checkIfExists: true)
                        .splitCsv(header: true)
                        .map { row ->
                            def meta = [
                                id         : row.sample,
                                input_type : row.input_type
                        ]
                        def h5ad_file = file(row.h5ad)
                        tuple(meta, h5ad_file)
                    }
                :
                    H5AD_CONVERSION.out.h5ad_filtered
            )
        CONVERT_MUDATA(
            ch_h5ad_selected,
            ch_vdj,
            ch_metadata_demuxafy,
            ch_metadata
        )
        ch_versions = ch_versions.mix(CONVERT_MUDATA.out.versions)
        ch_mudata = CONVERT_MUDATA.out.h5mu
    } else {
        ch_mudata = channel.empty()
    }

    //
    // SUBWORKFLOW: Run quality filtering on the concatenated h5ad files
    //
    // Da togliere questa cosa ch_rds_selected, se counts, canale vuoto tanto non faro' la parte dei doppietti
    def ch_rds_selected = params.counts ? H5AD_CONVERSION.out.rds_cellbender : H5AD_CONVERSION.out.rds_concat
    if ( !params.skip_qcfilters ) {
        DOUBLETS_QUALITYFILTERING (
            ch_rds_selected,
            CONVERT_MUDATA.out.h5mu,
            params.mt_threshold,
            params.min_umi_gex,
            params.max_umi_gex,
            params.min_genes_gex,
            params.max_genes_gex,
            params.min_cells_gex,
            params.min_features_adt,
            params.min_counts_adt
        )
        ch_versions = ch_versions.mix(DOUBLETS_QUALITYFILTERING.out.ch_versions)
        ch_h5mu_filtered = DOUBLETS_QUALITYFILTERING.out.h5mu
    } else {
        ch_h5mu_filtered = CONVERT_MUDATA.out.h5mu
    }

    //
    // SUBWORKFLOW: Run normalization on the concatenated h5ad files
    //
    ch_cellcycle_file = params.cell_cycle_file ?
        file(params.cell_cycle_file, checkIfExists: true) :
        channel.empty()

    // Make raw h5ad optional for reclustering workflows
    ch_h5ad_raw = params.h5ad_matrix ?
        channel.fromPath("${projectDir}/assets/EMPTY").map { file -> [[:], file] } :
        H5AD_CONVERSION.out.h5ad_raw

    NORMALIZATION_AND_HVG (
        ch_h5mu_filtered,
        ch_h5ad_raw,
        ch_cellcycle_file,
        params.n_pcs,
        params.n_neighbors,
        params.min_dist
    )
    ch_versions = ch_versions.mix(NORMALIZATION_AND_HVG.out.ch_versions)

    //
    // SUBWORKFLOW: Run cell annotation on the concatenated h5ad files
    //
    ch_input_model = params.input_model ? file(params.input_model, checkIfExists: true) : channel.empty()

    if ( params.input_model ) {
        CELL_ANNOTATION (
            NORMALIZATION_AND_HVG.out.h5mu,
            ch_input_model
        )
        ch_versions = ch_versions.mix(CELL_ANNOTATION.out.versions)
        ch_mu5ad = CELL_ANNOTATION.out.h5mu
        cell_annotation_meta_ch = CELL_ANNOTATION.out.metadata
    } else {
        ch_mu5ad = NORMALIZATION_AND_HVG.out.h5mu
        cell_annotation_meta_ch = channel.empty()
    }

    //
    // SUBWORKFLOW: Run ATAC preprocessing
    //
    atac_out_h5ad = channel.empty()

    if (params.aligner == "cellrangerarc") {
        blacklist_path = params.blacklist_path ? \
                         channel.value(file(params.blacklist_path, checkIfExists: true)) : \
                         channel.empty()

        ATAC_PREPROCESSING (
            ch_transformed_fragments_channel,
            ch_transformed_fragments_index_channel,
            params.tss_threshold,
            params.min_fragments_counts,
            params.max_fragments_counts,
            params.n_features_atac,
            params.frac_dup,
            params.peaks_frac,
            params.n_comps_atac,
            params.n_neighbors_atac,
            params.n_clusters_atac,
            blacklist_path,
            cell_annotation_meta_ch
        )
        atac_out_h5ad = ATAC_PREPROCESSING.out.h5ad
        ch_versions = ch_versions.mix(ATAC_PREPROCESSING.out.ch_versions)
    }

    //
    // SUBWORKFLOW: Run integration for GEX and ADT indipendently and jointly
    //

    INTEGRATION_MODALITIES (
        ch_mu5ad,
        atac_out_h5ad,
        params.n_neighbors_harmony,
        params.min_dist_harmony,
        params.integration_var
    )
    ch_versions = ch_versions.mix(INTEGRATION_MODALITIES.out.ch_versions)

    //
    // MODULES: Run clustering for GEX
    //
    CLUSTERING (
        INTEGRATION_MODALITIES.out.h5mu_out,
        params.resolution_min,
        params.resolution_max,
        params.top_n_markers
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

    //
    // MODULES: Plot clustree graph
    //
    CLUSTREE (
        CLUSTERING.out.metadata_final
    )
    ch_versions = ch_versions.mix(CLUSTREE.out.versions)

    // Handling multiple resolutions
    if ( params.resolution ) {
        resolution_ch = channel.fromList(params.resolution.toString().split(',').flatten())

        //
        // MODULES: Enrichment on marker genes for a selected resolution
        //
        if ( params.enrich_collection ){
            ch_enrich_collection = channel.fromList(params.enrich_collection.split(',').flatten())
            resolution_ch
                .combine( ch_enrich_collection )
                .map{ res, coll -> [["res": res, "coll": coll], res, coll] }
                .set { ch_res_enrich }

            ENRICH_MARKERS (
                CLUSTERING.out.ranked_genes.collect(),
                ch_res_enrich
            )
            ch_versions = ch_versions.mix(ENRICH_MARKERS.out.versions)
        }
    }

    //
    // MODULES: Plot custom genelist
    //
    if ( params.custom_geneset ) {
        ch_custom_geneset = channel.fromList(params.custom_geneset.split(',').flatten())

        if ( params.resolution ) {
            resolution_ch
                .combine( ch_custom_geneset )
                .map{ res, genes -> [["res": res, "genes": genes], res, genes] }
                .set { ch_res_geneset }
        } else {
            // if no resolution is provided, use 100 as fake resolution
            fake_res = 100
            ch_res_geneset = ch_custom_geneset.map { genes ->
                [["res": fake_res, "genes": genes], fake_res, genes]
            }
        }
        CUSTOM_GENES (
            CLUSTERING.out.h5mu.collect(),
            ch_res_geneset
        )
        ch_versions = ch_versions.mix(CUSTOM_GENES.out.versions)
    }

    '''
    //
    // MODULE: Run differential analysis
    //
    DIFFERENTIAL_ANALYSIS (
        CLUSTERING.out.h5mu
    )
    ch_versions = ch_versions.mix(DIFFERENTIAL_ANALYSIS.out.versions)
    '''

    if (params.resolution) {

        resolution_ch = channel.fromList(params.resolution.toString().split(',').flatten())

        DIFFERENTIAL_ABUNDANCE(
            CLUSTERING.out.h5mu
                .combine(ch_diff_abundance_comparisons)
                .combine(resolution_ch)
        )
        if (DIFFERENTIAL_ABUNDANCE.out.versions) {
            ch_versions = ch_versions.mix(DIFFERENTIAL_ABUNDANCE.out.versions)
        }

    }

    if ( params.resolution ) {

        ch_resolution = channel.fromList(params.resolution.toString().split(',').flatten())

        PSEUDOBULK_ANALYSIS(
            CLUSTERING.out.h5mu,
            ch_resolution,
            ch_pseudobulk_group,
            ch_pseudobulk_comparisons,
            ch_pseudobulk_formula,
            ch_pseudobulk_fdr
        )
        if (PSEUDOBULK_ANALYSIS.out.versions) {
            ch_versions = ch_versions.mix(PSEUDOBULK_ANALYSIS.out.versions)
        }
    }

    // Cell to cell interaction
    if ( params.resolution ) {

        ch_resolution = channel.fromList(params.resolution.toString().split(',').flatten())

        CLUSTERING.out.h5mu
            .combine(ch_liana_method)
            .combine(ch_liana_resource)
            .combine(ch_resolution)
            .set { cell_interaction_input }

        CELL_INTERACTION(
            cell_interaction_input
        )

        if (CELL_INTERACTION.out.versions) {
            ch_versions = ch_versions.mix(CELL_INTERACTION.out.versions)
        }

    }


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'scrnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    if (!params.skip_multiqc) {
        //
        // MODULE: MultiQC
        //
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        def ch_multiqc_custom_methods_description = multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'scrnaseq'],
                    files,
                    multiqc_config
                        ? file(multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    } else {
        ch_multiqc_report = channel.empty()
    }

    emit:
    multiqc_report = ch_multiqc_report           // channel: [ path(multiqc_report.html) ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
