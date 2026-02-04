process CLUSTERING  {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/nfdata/sc_rnaseq:v1.0.1'

    input:
    tuple val(meta), path(input_h5mu)
    val resolution_min
    val resolution_max
    val top_n_markers

    output:
    tuple val(meta), path("matrix.clustered*.h5mu")   , emit: h5mu
    tuple val(meta), path("final_metadata.csv")       , emit: metadata_final
    path "cluster_id*.pdf"                            , emit: clusters
    path "ranked_genes.xlsx"                          , emit: ranked_genes
    path "top_*_markers_heatmap.pdf"                  , emit: heatmaps
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    clustering.py -ad $input_h5mu -min_res $resolution_min -max_res $resolution_max -n $top_n_markers

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        clustering.py: \$(clustering.py --version 2> /dev/null | grep -v scanpy)
    END_VERSIONS
    """

    stub:
    """
    touch matrix_clustered.h5mu
    touch final_metadata.csv
    touch cluster_id_all.pdf

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        clustering.py: \$(clustering.py --version 2> /dev/null | grep -v scanpy)
    END_VERSIONS
    """


}
