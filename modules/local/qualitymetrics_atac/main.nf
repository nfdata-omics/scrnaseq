process QUALITY_FILTERING_ATAC  {
    tag "$meta.id"
    label 'process_high_memory'

    container 'docker.io/nfdata/snapatac:v1.0.0'

    input:
    tuple val(meta), path(input_fragment_file,stageAs: "?/*")
    tuple val(meta2), path(input_fragment_index_file,stageAs: "?/*")
    val tss_threshold
    val min_fragments_counts
    val max_fragments_counts
    path blacklist_path
    path genome_annotation

    output:
    tuple val(meta), path("*.filtered_atac.h5ad"), emit: h5ad
    path "FragSizeDist_all_samples.pdf", emit: fragment_size_distribution, optional: true
    path "QC_Histograms_all_samples.pdf", emit: qc_histograms, optional: true
    path "TSSE_score_sample_*", emit: TSS_score, optional: true
    path "TSSE_score_all_samples.pdf", emit: tss_signal, optional: true
    path "cell_counts.csv", emit: cell_counts, optional: true
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fragment_files = input_fragment_file instanceof List ? input_fragment_file : [input_fragment_file]
    def fragment_index_files = input_fragment_index_file instanceof List ? input_fragment_index_file : [input_fragment_index_file]
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export XDG_CACHE_HOME=/tmp

    qualitymetricsfilters_atac.py  -fr ${fragment_files.join(' ')} -fri ${fragment_index_files.join(' ')}  -id ${meta.collect{ it.id }.join(' ')} -t $tss_threshold -mif $min_fragments_counts -maf $max_fragments_counts -b $blacklist_path -g $genome_annotation
    
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        qualitymetricsfilters_atac.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix.filtered_atac.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualitymetricsfilters_atac.py --version >> versions.yml
    END_VERSIONS
    """
}
