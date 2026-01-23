process QUALITY_FILTERING_ATAC  {
    tag "$meta.id"
    label 'process_high_memory'

    container 'docker.io/nfdata/snapatac:v1.0.0'


    input:
    tuple val(meta), path (input_fragment_file,stageAs: "?/*")
    tuple val(meta), path (input_fragment_index_file,stageAs: "?/*")
    val tss_threshold
    val min_fragments_counts
    val max_fragments_counts
    path blacklist_path

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
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export XDG_CACHE_HOME=/tmp
    export POOCH_CACHE_DIR=$PWD/.pooch

    python /home/camilla.callierotti/NFG.2118962_NeuroMuscolar_Carra_CC/scrnaseq/modules/local/qualitymetrics_atac/resources/usr/bin/qualitymetricsfilters_atac.py \
    -fr ${input_fragment_file instanceof List ? input_fragment_file.join(' ') : input_fragment_file} \
    -fri ${input_fragment_index_file instanceof List ? input_fragment_index_file.join(' ') : input_fragment_index_file} \
    -id ${meta instanceof List ? meta.collect{ it.id }.join(' ') : meta.id} \
    -t $tss_threshold -mif $min_fragments_counts -maf $max_fragments_counts -b $blacklist_path
    
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
