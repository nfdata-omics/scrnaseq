process ENRICH_MARKERS  {
    tag "${meta.res}_${meta.coll}"
    label 'process_single'

    container 'docker.io/nfdata/clusterprofiler:v4.14.4'

    input:
    path ranked_genes
    tuple val(meta), val(resolution), path(collection)

    output:
    path "enrich_Leiden_*.xlsx", emit: enriched_markers
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp


    enrich_markergenes.R -b $ranked_genes $resolution $collection

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        enrich_markergenes.R --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch enrich_Leiden_.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        enrich_markergenes.R --version >> versions.yml
    END_VERSIONS
    """

}
