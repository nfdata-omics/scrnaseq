process NORMALIZATION   {
    tag "$meta.id"
    label 'process_single'

    container = 'quay.io/biocontainers/scirpy:0.20.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(input_h5mu)

    output:
    tuple val(meta), path("*.norm.h5mu"), emit: h5mu
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    normalization.py -ad $input_h5mu

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
    END_VERSIONS
    normalization.py --version >> versions.yml
    """

    stub:
    """
    touch matrix.norm.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    normalization.py --version >> versions.yml
    """

}