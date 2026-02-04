process CONVERT_MUDATA  {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/scirpy:0.20.1--pyhdfd78af_0'

    input:
    tuple val(meta),  path(input_h5ad)
    tuple val(meta2), path(input_vdj)
    tuple val(meta3), path(demultiplexing_doublets)
    path (metadata_file)

    output:
    tuple val(meta), path("*.mudata.h5mu") , emit: h5mu
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ai = input_vdj ? "-ai $input_vdj" : ''
    def csv = demultiplexing_doublets ? "-csv $demultiplexing_doublets" : ''
    def meta_file = metadata_file ? "-meta $metadata_file" : ''
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    convert.py -ad $input_h5ad $ai $csv $meta_file

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        convert.py: \$(convert.py --version 2> /dev/null | grep -v scanpy)
    END_VERSIONS
    """

    stub:
    """
    touch matrix.mudata.h5mu

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        convert.py: \$(convert.py --version 2> /dev/null | grep -v scanpy)
    END_VERSIONS
    """
}
