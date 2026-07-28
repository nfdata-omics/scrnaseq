//
// Subworkflow with functionality specific to the nfdata-omics/scrnaseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    _input            //  string: Path to input samplesheet
    h5ad_matrix       //  string: Path to input h5ad matrix file
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/scrnaseq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/scrnaseq/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    if (h5ad_matrix) {
        ch_h5ad_matrix = h5ad_matrix ? channel.fromPath( h5ad_matrix, checkIfExists: true ) : channel.empty()
        ch_samplesheet = channel.empty()
    }
    else {
        //
        // Create channel from input file provided through params.input
        //
        PARSE_SAMPLESHEET(params.input)
        ch_samplesheet = PARSE_SAMPLESHEET.out.samplesheet
        ch_h5ad_matrix = channel.empty()
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
    h5ad_matrix = ch_h5ad_matrix
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PARSE THE INPUT SAMPLESHEET
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PARSE_SAMPLESHEET {

    take:
    samplesheet // path: Input samplesheet

    main:
    if (!samplesheet) {
        error("No input samplesheet was provided. Please specify `--input`.")
    }
    def samplesheet_rows = samplesheetToList(samplesheet, "${projectDir}/assets/schema_input.json")
    validateInputSamplesheet(samplesheet_rows)

    ch_samplesheet_rows = channel.fromList(samplesheet_rows)

    ch_fastq = ch_samplesheet_rows
        .filter { _meta, fastq_1, _fastq_2, _processed_data, _unfiltered_data ->
            hasSamplesheetValue(fastq_1)
        }
        .map { meta, fastq_1, fastq_2, _processed_data, _unfiltered_data ->
            def reads = [fastq_1]
            if (hasSamplesheetValue(fastq_2)) {
                reads << fastq_2
            }
            if (meta.feature_type == 'atac') {
                reads << file(meta.fastq_barcode, checkIfExists: true)
            }

            [
                meta.id,
                meta.feature_type,
                meta + [
                    input_type: 'fastq',
                    single_end: !hasSamplesheetValue(fastq_2)
                ],
                reads
            ]
        }
        .groupTuple(by: [0, 1])
        .map { id, _feature_type, metas, fastqs ->
            def (meta, validated_fastqs) = validateFastqSample([id, metas, fastqs])
            [meta, validated_fastqs.flatten()]
        }

    ch_preprocessed = ch_samplesheet_rows
        .filter { _meta, _fastq_1, _fastq_2, processed_data, _unfiltered_data ->
            hasSamplesheetValue(processed_data)
        }
        .flatMap { meta, _fastq_1, _fastq_2, processed_data, unfiltered_data ->
            def entries = [
                [meta + [input_type: 'filtered'], [processed_data]]
            ]
            if (hasSamplesheetValue(unfiltered_data)) {
                entries << [meta + [input_type: 'raw'], [unfiltered_data]]
            }
            entries
        }

    ch_fastq
        .mix(ch_preprocessed)
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Retrieve the aligner-specific protocol based on the specified protocol.
// Returns a map containing the protocol and any optional aligner configuration.
def getProtocol(workflow, log, aligner, protocol) {
    def jsonSlurper = new groovy.json.JsonSlurper()
    def protocols = jsonSlurper.parseText(new File("${workflow.projectDir}/assets/protocols.json").text)
    def protocol_aligner = aligner == 'cellrangermulti' ? 'cellranger' : aligner
    def aligner_map = protocols[protocol_aligner]

    if (aligner_map.containsKey(protocol)) {
        return aligner_map[protocol]
    }

    log.warn("Protocol '${protocol}' not recognized by the pipeline. Passing on the protocol to the aligner unmodified.")
    return [protocol: protocol]
}
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    // Validate cellranger_multi_barcodes if provided and aligner is cellrangermulti
    if (params.aligner == 'cellrangermulti' && params.cellranger_multi_barcodes  && params.input) {
        validateCellrangerMultiBarcodes()
    }
}

//
// Validate cellranger_multi_barcodes samplesheet for uniqueness and conditional requirements
//
def validateCellrangerMultiBarcodes() {
    def cellranger_multi_barcodes = file(params.cellranger_multi_barcodes).splitCsv(header: true)

    // Get unique samples from input samplesheet for cross-validation
    def inputSamples = file(params.input).splitCsv(header: true).collect { row -> row.sample }.toSet()

    // Check that at least one barcode column is provided for each row
    // and that each sample uses only one type of barcode
    def rowsWithoutBarcodes = []
    def sampleBarcodeTypes = [:]
    def barcodeSamples = [] as Set
    cellranger_multi_barcodes.eachWithIndex { row, idx ->
        def multiplexed_sample_id = row.multiplexed_sample_id
        def rowNum = idx + 2 // +2 for 1-based indexing and header row

        // Collect unique sample names for cross-validation
        barcodeSamples << row.sample

        def barcodeTypes = []
        if (row.probe_barcode_ids) barcodeTypes << 'probe_barcode_ids'
        if (row.cmo_ids)           barcodeTypes << 'cmo_ids'
        if (row.ocm_ids)           barcodeTypes << 'ocm_ids'

        if (barcodeTypes.isEmpty()) {
            rowsWithoutBarcodes << [row: rowNum, multiplexed_sample_id: multiplexed_sample_id]
        }
        sampleBarcodeTypes[multiplexed_sample_id] = [types: barcodeTypes.toSet(), row: rowNum]
    }

    // Validate that at least one barcode identifier is populated in each row
    if (rowsWithoutBarcodes) {
        def errorDetails = rowsWithoutBarcodes.collect { missing -> "row ${missing.row} (${missing.multiplexed_sample_id})" }.join(', ')
        error("Please check cellranger_multi_barcodes samplesheet -> " +
              "The following rows have no barcode identifiers: ${errorDetails}. " +
              "Each row must have exactly one of: 'probe_barcode_ids', 'cmo_ids', or 'ocm_ids'.")
    }

    // Validate that no more than one barcode identifier is populated in each row
    def samplesWithMixedBarcodes = sampleBarcodeTypes.findAll { _multiplexed_sample_id, info -> info.types.size() > 1 }
    if (samplesWithMixedBarcodes) {
        def errorMsg = samplesWithMixedBarcodes.collect { multiplexed_sample_id, info ->
            "'${multiplexed_sample_id}' (row ${info.row}) uses multiple barcode types: ${info.types.join(', ')}"
        }.join('; ')
        error("Please check cellranger_multi_barcodes samplesheet -> " +
              "Each multiplexed_sample_id should use only one type of barcode identifier. ${errorMsg}")
    }

    // Validate that samples in cellranger_multi_barcodes exist in the input samplesheet
    def unknownSamples = barcodeSamples - inputSamples
    if (unknownSamples) {
        error("Please check cellranger_multi_barcodes samplesheet -> " +
              "The following sample(s) do not exist in the input samplesheet: ${unknownSamples.join(', ')}. " +
              "The 'sample' column in cellranger_multi_barcodes must match 'sample' values in the input samplesheet.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    // samplesheetToList returns one list per samplesheet row, with metadata in
    // the first position. Validate constraints that span more than one row
    // before constructing the input channel.
    def rows_by_sample = [:].withDefault { [] }

    input.eachWithIndex { row, index ->
        def (meta, fastq_1, fastq_2, processed_data, unfiltered_data) = row
        def sample = meta.id
        def row_number = index + 2

        if (!hasSamplesheetValue(meta.feature_type)) {
            meta.feature_type = 'gex'
        }

        def valid_feature_types = params.aligner == 'cellrangerarc'
            ? ['gex', 'atac']
            : params.aligner == 'cellrangermulti'
                ? ['gex', 'vdj', 'ab', 'crispr', 'cmo']
                : ['gex']
        if (meta.feature_type !in valid_feature_types) {
            error("Please check input samplesheet -> feature_type '${meta.feature_type}' is not supported by aligner '${params.aligner}' for sample '${sample}'. Allowed values: ${valid_feature_types.join(', ')}.")
        }
        if (params.aligner == 'cellrangerarc' && meta.feature_type == 'atac' && !hasSamplesheetValue(meta.fastq_barcode)) {
            error("Please check input samplesheet -> Cell Ranger ARC requires a barcode FASTQ for ATAC input: ${sample}.")
        }

        rows_by_sample[sample] << [
            row_number     : row_number,
            has_fastq      : hasSamplesheetValue(fastq_1) || hasSamplesheetValue(fastq_2),
            processed_data : processed_data,
            unfiltered_data: unfiltered_data
        ]

        if (hasSamplesheetValue(processed_data)) {
            validatePreprocessedInput(processed_data, sample, 'processed_data')
        }
        if (hasSamplesheetValue(unfiltered_data)) {
            validatePreprocessedInput(unfiltered_data, sample, 'unfiltered_data')
        }
    }

    rows_by_sample.each { sample, rows ->
        def fastq_rows = rows.findAll { row -> row.has_fastq }
        def processed_rows = rows.findAll { row -> hasSamplesheetValue(row.processed_data) }
        def unfiltered_rows = rows.findAll { row -> hasSamplesheetValue(row.unfiltered_data) }

        if (fastq_rows && (processed_rows || unfiltered_rows)) {
            error("Please check input samplesheet -> Sample '${sample}' provides both FASTQ and preprocessed input. Choose only one input form.")
        }
        if (unfiltered_rows && !processed_rows) {
            error("Please check input samplesheet -> Sample '${sample}' provides unfiltered_data without the required processed_data.")
        }
        if (processed_rows.size() > 1) {
            def row_numbers = processed_rows.collect { row -> row.row_number }.join(', ')
            error("Please check input samplesheet -> Sample '${sample}' provides processed_data more than once (rows ${row_numbers}).")
        }
        if (unfiltered_rows.size() > 1) {
            def row_numbers = unfiltered_rows.collect { row -> row.row_number }.join(', ')
            error("Please check input samplesheet -> Sample '${sample}' provides unfiltered_data more than once (rows ${row_numbers}).")
        }
    }

    return input
}

//
// Validate grouped FASTQ rows for one sample
//
def validateFastqSample(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Check whether an optional samplesheet value is populated
//
def hasSamplesheetValue(value) {
    return value != null &&
        (!(value instanceof Collection) || !value.isEmpty()) &&
        value.toString().trim()
}

//
// Validate supported preprocessed inputs: a Cell Ranger HDF5 file or a
// MEX/MTX directory containing the three required matrix components.
//
def validatePreprocessedInput(input, sample, column) {
    def input_path = input as java.nio.file.Path

    if (java.nio.file.Files.isRegularFile(input_path)) {
        if (!input_path.fileName.toString().toLowerCase().endsWith('.h5')) {
            error("Please check input samplesheet -> ${column} for sample '${sample}' must be an .h5 file or a MEX/MTX directory: ${input_path}")
        }
        if (!java.nio.file.Files.isReadable(input_path)) {
            error("Please check input samplesheet -> ${column} for sample '${sample}' is not readable: ${input_path}")
        }
        return
    }

    if (java.nio.file.Files.isDirectory(input_path)) {
        if (!java.nio.file.Files.isReadable(input_path)) {
            error("Please check input samplesheet -> ${column} directory for sample '${sample}' is not readable: ${input_path}")
        }

        def required_files = [
            matrix  : ['matrix.mtx', 'matrix.mtx.gz'],
            barcodes: ['barcodes.tsv', 'barcodes.tsv.gz'],
            features: ['features.tsv', 'features.tsv.gz']
        ]
        def missing_files = required_files.findAll { name, alternatives ->
            !alternatives.any { filename ->
                def candidate = input_path.resolve(filename)
                java.nio.file.Files.isRegularFile(candidate) && java.nio.file.Files.isReadable(candidate)
            }
        }.keySet()

        if (missing_files) {
            error("Please check input samplesheet -> ${column} for sample '${sample}' is not a valid MEX/MTX directory. Missing readable component(s): ${missing_files.join(', ')}.")
        }
        return
    }

    error("Please check input samplesheet -> ${column} for sample '${sample}' must be an .h5 file or a MEX/MTX directory: ${input_path}")
}
//
// cellrangerarc structure for samplesheet channel
//
def cellrangerarcStructure(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Validate that the property "feature_type" has a valid value for Cell Ranger ARC
    def valid_feature_types = ["gex", "atac"]
    def feature_type_ok = metas.collect { meta -> meta.feature_type }.unique().every { feature_type -> feature_type in valid_feature_types }
    if (!feature_type_ok) {
        error("Please check input samplesheet -> For cellrangerarc, 'feature_type' can only be 'gex' or 'atac'.")
    }

    // Define a new common meta for all the fastqs in this channel instance
    def sampleMeta = metas[0].clone()
    sampleMeta.remove("feature_type")

    // Create a list with all the feature types expected by Cell Ranger ARC
    def sampletypes = metas.collect { meta -> meta.feature_type }

    // Create a list with all the base name of the fastq files
    def subsamples = fastqs.collect { fastq ->
        def match = (fastq[0].baseName =~ /^(.*?)_S\d+_L\d+_R\d+_\d+\.fastq(\.gz)?$/)
        if (!match) {
            error("Filename does not follow the expected FASTQ filename convention (SampleName_S1_L001_R1_001.fastq.gz): ${fastq[0]}")
        }
        return match[0][1]
    }

    return [ sampleMeta, sampletypes, subsamples, fastqs.flatten() ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        } else {
            return null
        }
    } else {
        return null
    }
}

//
// iGenomes GTF annotations with spaces in the GTF source column (e.g. NCBI RefSeq "Curated Genomic")
// are incompatible with Cell Ranger 10 mkref; opt-in per genome via gtf_source_has_spaces.
//
def gtfSourceFixNeeded(aligner, genome, genomes, gtf) {
    def genome_entry = genomes && genome ? genomes[genome] : null
    def cellranger_aligner = aligner in ['cellranger', 'cellrangerarc', 'cellrangermulti']
    def gtf_flagged = genome_entry?.gtf_source_has_spaces as Boolean
    def gtf_from_genome = gtf == genome_entry?.gtf
    return cellranger_aligner && gtf_flagged && gtf_from_genome
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
