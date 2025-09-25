process CLUSTREE  {
    tag "$meta.id"
    label 'process_single'

    container = 'docker.io/nfdata/sc-rnaseq-r:v1.0.0'

    input:
    tuple val(meta), path(cluster_id)

    output:
    path "clustree_plot.png", emit: clustree
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp


    clustree.R $cluster_id

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        clustree.R --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch clustree_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustree.R --version >> versions.yml
    END_VERSIONS
    """

}
