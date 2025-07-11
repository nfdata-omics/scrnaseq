process MOFA_INTEGRATION {
    tag "$meta.id"
    label 'process_medium'

    container = 'docker.io/nfdata/muon-sc_rnaseq:v1.0.2'

    input:
    tuple val(meta), path(input_h5mu)

    output:
    tuple val(meta), path("*.mofa.h5mu"), emit: h5mu
    path "umap_mofa.png", emit: graph_UMAP_mofa
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    mofa_integration.py -ad $input_h5mu

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        mofa_integration.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix.mofa.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mofa_integration.py --version >> versions.yml
    END_VERSIONS
    """

}
