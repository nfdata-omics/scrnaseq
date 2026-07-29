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
include { ALIGNMENT                                         } from '../subworkflows/local/alignment'
include { FASTQC_CHECK                                      } from '../subworkflows/local/fastqc'
include { MTX_TO_H5AD                                       } from '../modules/local/mtx_to_h5ad'
include { H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA } from '../subworkflows/nf-core/h5ad_removebackground_barcodes_cellbender_anndata'
include { H5AD_CONVERSION                                   } from '../subworkflows/local/h5ad_conversion'
include { ATAC_PREPROCESSING                                } from '../subworkflows/local/preprocessing_atac'
include { CONCATENATE_VDJ                                   } from '../modules/local/concatenate_vdj'
include { QUALITY_METRICS_VDJ                               } from '../modules/local/qualitymetrics_vdj'
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
    ch_samplesheet              // channel: [ meta, files ] from samplesheet
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

    ch_fastq = ch_samplesheet.filter { meta, _files -> meta.input_type == 'fastq' }

    protocol_config = getProtocol(workflow, log, params.aligner, params.protocol)
    if (protocol_config['protocol'] == 'auto' && params.aligner !in ["cellranger", "cellrangerarc", "cellrangermulti"]) {
        error "Only cellranger supports `protocol = 'auto'`. Please specify the protocol manually!"
    }

    // general input and params
    ch_motifs               = motifs           ? file(motifs, checkIfExists: true)           : []
    // Warn if both GTF and GFF files are provided
    if (gtf && gff) {
        log.warn("Both GTF and GFF files are provided. GTF file will be used.")
    }

    // samplesheet - this is passed to the MTX conversion functions to add metadata to the
    // AnnData objects.
    ch_input = params.input                ? file(params.input, checkIfExists: true)    : []
    ch_h5ad_matrix = params.h5ad_matrix    ? file(params.h5ad_matrix, checkIfExists: true): []

    //cellranger params
    ch_cellranger_index = cellranger_index ? file(cellranger_index, checkIfExists: true) : []

    //cellrangermulti params
    cellranger_vdj_index = cellranger_vdj_index             ? file(cellranger_vdj_index, checkIfExists: true)             : []
    ch_multi_samplesheet = params.cellranger_multi_barcodes ? file(params.cellranger_multi_barcodes, checkIfExists: true) : []

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

    //
    // Run alignment pipeline
    //
    ALIGNMENT(
        ch_samplesheet,
        ch_genome_fasta,
        ch_genome_gtf,
        ch_cellranger_index,
        protocol_config['protocol'],
        ch_motifs,
        ch_cellrangerarc_config,
        cellranger_vdj_index,
        ch_multi_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.multiqc_files.flatten())

    //
    // MODULE: Convert mtx matrices to h5ad
    //
    def ch_matrix_inputs = ALIGNMENT.out.mtx_matrices.map { meta, matrix_paths ->
        tuple(meta, selectMatrixInput(matrix_paths, meta))
    }
    MTX_TO_H5AD (
        ch_matrix_inputs
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
        ch_input
    )
    ch_versions = ch_versions.mix(H5AD_CONVERSION.out.ch_versions)

    //
    // MODULE: Concat vdj samples and save as h5ad format
    //

    def vdj_file = channel.empty()
    if (params.vdj) {

        // Case 1: VDJ provided as external input
        ch_vdj_input = channel
            .fromPath(params.vdj, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample]
                def contig_file = file(row.csv, checkIfExists: true)

                tuple(meta, contig_file)
            }
            .collect(flat: false)
            .map { rows ->
                tuple(
                    rows.collect { it[0] },
                    rows.collect { it[1] }
                )
            }

        CONCATENATE_VDJ(ch_vdj_input)

        ch_versions = ch_versions.mix(
            CONCATENATE_VDJ.out.versions
        )

        vdj_file = CONCATENATE_VDJ.out.h5ad

    } else if (
        params.aligner == "cellrangermulti" &&
        params.include_vdj
    ) {

        // Case 2: VDJ generated by Cell Ranger Multi
        CONCATENATE_VDJ(
            ALIGNMENT.out.vdj_file
        )

        ch_versions = ch_versions.mix(
            CONCATENATE_VDJ.out.versions
        )

        vdj_file = CONCATENATE_VDJ.out.h5ad

    }
    // If no VDJ input is available, send a dummy tuple to CONVERT_MUDATA
    vdj_file = vdj_file.ifEmpty {
        tuple([id: 'dummy'], [])
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
        def ch_h5ad_selected = params.h5ad_matrix ?
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
        CONVERT_MUDATA(
            ch_h5ad_selected,
            vdj_file,
            ch_metadata_demuxafy,
            ch_metadata
        )
        ch_versions = ch_versions.mix(CONVERT_MUDATA.out.versions)
        ch_mudata = CONVERT_MUDATA.out.h5mu
    } else {
        ch_mudata = channel.empty()
    }

    //
    // MODULE: Run quality filtering on the vdj concatenated h5ad files
    //

    if (params.vdj || params.include_vdj) {

        QUALITY_METRICS_VDJ(
            CONVERT_MUDATA.out.h5mu
        )

        ch_versions = ch_versions.mix(
            QUALITY_METRICS_VDJ.out.versions
        )
    }

    //
    // SUBWORKFLOW: Run quality filtering on the concatenated h5ad files
    //
    if ( !params.skip_qcfilters ) {
        DOUBLETS_QUALITYFILTERING (
            H5AD_CONVERSION.out.rds_concat,
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
        genome_annotation_path = params.genome_annotation_path ? \
                         channel.value(file(params.genome_annotation_path, checkIfExists: true)) : \
                         channel.empty()


        ATAC_PREPROCESSING (
            ALIGNMENT.out.fragments_file,
            ALIGNMENT.out.fragments_index,
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
            genome_annotation_path,
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
        params.skip_harmony,
        params.integration_var,
        params.skip_integration
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

/*
 * Cell Ranger modules currently expose all files below their output directory,
 * whereas samplesheet restart entries contain one logical matrix artifact.
 * Prefer the H5 representation when both H5 and MEX are present in the legacy
 * aligner output. A supplied MEX directory is used directly.
 */
def selectMatrixInput(matrix_paths, meta) {
    def paths = matrix_paths instanceof Collection ? matrix_paths as List : [matrix_paths]
    def h5_paths = paths.findAll { path ->
        path.name.toLowerCase().endsWith('.h5')
    }

    if (h5_paths.size() == 1) {
        return h5_paths.first()
    }
    if (h5_paths.size() > 1) {
        error("MTX_TO_H5AD received multiple H5 matrices for sample '${meta.id}' (${meta.input_type}): ${h5_paths*.name.join(', ')}")
    }

    def mex_paths = paths.findAll { path ->
        isMexDirectory(path)
    }
    if (mex_paths.size() == 1) {
        return mex_paths.first()
    }
    if (mex_paths.size() > 1) {
        error("MTX_TO_H5AD received multiple MEX directories for sample '${meta.id}' (${meta.input_type}): ${mex_paths*.name.join(', ')}")
    }

    error("MTX_TO_H5AD requires one H5 matrix or one complete MEX directory for sample '${meta.id}' (${meta.input_type}).")
}

def isMexDirectory(path) {
    def input_path = path as java.nio.file.Path
    if (!java.nio.file.Files.isDirectory(input_path)) {
        return false
    }

    def required_components = [
        ['matrix.mtx', 'matrix.mtx.gz'],
        ['barcodes.tsv', 'barcodes.tsv.gz'],
        ['features.tsv', 'features.tsv.gz']
    ]
    return required_components.every { alternatives ->
        alternatives.any { filename ->
            java.nio.file.Files.isRegularFile(input_path.resolve(filename))
        }
    }
}
