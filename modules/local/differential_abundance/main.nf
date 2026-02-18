process DIFFERENTIAL_ABUNDANCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container 'docker.io/nfdata/scverse-pertpy:v1.0.0-scanpy1.12'

    input:
    tuple val(meta), path(h5mu)
    val target
    val reference
    val column

    output:
    path "diff_abundance/*", emit: diff_abund_results

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    
    python3 ${moduleDir}/resources/usr/bin/differential_abundance.py \
        --mdata $h5mu \
        --target $target \
        --reference $reference \
        --column_to_test $column
    """

}
