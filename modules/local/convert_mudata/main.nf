process CONVERT_MUDATA  {
    tag "$meta.id"
    label 'process_medium'

    container = 'quay.io/biocontainers/scirpy:0.20.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(input_h5ad)
    tuple val(meta), path(input_vdj)
    tuple val(meta), path (demultiplexing_doublets)
    //tuple val(meta), path(input_h5ad_atac)

    output:
    tuple val(meta), path("*.mudata.h5mu") , emit: h5mu
    path "versions.yml",  emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def ai = input_vdj ? "-ai $input_vdj" : ''
    def csv = demultiplexing_doublets ? "-csv $demultiplexing_doublets" : ''

    //convert.py -ad $input_h5ad $ai -at $input_h5ad_atac

    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    convert.py -ad $input_h5ad $ai $csv

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        convert.py --version >> versions.yml
    END_VERSIONS   
    """

    stub:
    """
    touch matrix.mudata.h5mu
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        convert.py --version >> versions.yml
    END_VERSIONS
    """
}
