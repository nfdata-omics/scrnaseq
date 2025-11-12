process MOFA_INTEGRATION {
    tag "$meta.id"
    label 'process_high'

    container = 'docker.io/nfdata/muon-sc_rnaseq:v.1.0.5'
    
    input:
    tuple val(meta), path(input_h5mu)
    tuple val(meta), path(input_h5ad_atac)

    output:
    tuple val(meta), path("*.mofa.h5mu"), emit: h5mu
    path "umap_mofa.png", emit: umap_mofa, optional: true
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    mofa_integration.py -ad $input_h5mu -at $input_h5ad_atac

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
