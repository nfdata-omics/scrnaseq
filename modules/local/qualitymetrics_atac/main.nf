process QUALITY_FILTERING_ATAC  {
    tag "$meta.id"
    label 'process_high'

    container = 'quay.io/biocontainers/snapatac2:2.8.0--py311h284d45d_1'
    

    input:
    tuple val(meta), path (input_fragment_file, stageAs: "?/*")
    tuple val(meta), path (input_fragment_index_file, stageAs: "?/*")
    val nucleosome_threshold
    val tss_threshold
    
    output:
    tuple val(meta), path("*.filtered_atac.h5ad"), emit: h5ad
    path "Peaks_Distribution_*.png", emit: peaks_distribution, optional: true
    path "Nucleosome_signal_*.png", emit: nucleosome_signal, optional: true
    path "TSS_score_*.png", emit: tss_signal, optional: true
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    //qualitymetricsfilters.py -ad $input_h5mu -d $doublets_csv -fr ${input_fragment_file.join(' ')} -id ${meta.collect{ it.id }.join(' ')} -f $mt_threshold -n $nucleosome_threshold -t $tss_threshold
    
    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export XDG_CACHE_HOME=/tmp

    qualitymetricsfilters_atac.py  -fr ${input_fragment_file.join(' ')} -fri $input_fragment_index_file  -id ${meta.collect{ it.id }.join(' ')} -n $nucleosome_threshold -t $tss_threshold
    
    
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
