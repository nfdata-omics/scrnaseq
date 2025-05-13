process DIFFERENTIAL_ANALYSIS  {
    tag "$meta.id"
    label 'process_single'

    container = 'docker.io/nfdata/muon-sc_rnaseq:v1.0.4'

    input:
    tuple val(meta), path(input_h5mu)

    output:
    path "differential_genes.xlsx", emit: xlsx
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp


    differential_analysis.py -ad $input_h5mu

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        differential_analysis.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch differential_genes.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        differential_analysis.py --version >> versions.yml
    END_VERSIONS
    """

}
