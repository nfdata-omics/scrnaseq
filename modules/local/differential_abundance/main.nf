process DIFFERENTIAL_ABUNDANCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container = 'docker.io/nfdata/muon-sc_rnaseq:v1.0.4'

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
    python3 \${moduleDir}/resources/usr/bin/differential_abundance.py \
        --mdata $h5mu \
        --target $target \
        --reference $reference \
        --column_to_test $column
    """

}