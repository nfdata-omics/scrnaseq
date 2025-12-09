process CLUSTERING  {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/nfdata/sc_rnaseq:v1.0.1'

    input:
    tuple val(meta), path(input_h5mu)
    val resolution

    output:
    tuple val(meta), path("*.clustered.h5mu")   , emit: h5mu
    tuple val(meta), path("final_metadata.csv") , emit: metadata_final
    path "cluster_id*.pdf"                      , emit: clusters
    path "ranked_genes.xlsx"                    , emit: ranked_genes
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    clustering.py -ad $input_h5mu -res $resolution

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        clustering.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix_clustered.h5ad
    touch final_metadata.csv
    touch cluster_id.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustering.py --version >> versions.yml
    END_VERSIONS
    """


}
