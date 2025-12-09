process HIGHLY_VARIABLE_GENES  {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/scirpy:0.20.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(input_h5mu)
    val n_pcs
    val n_neighbors
    val min_dist

    output:
    tuple val(meta), path("*.hvg.h5mu"), emit: h5mu
    path "umap_coordinates.csv"        , emit: umap
    path "umap_plot_*.pdf"             , emit: graph_umap
    path "pca_GEX.pdf"                 , emit: pca_gex
    path "pca_elbow.pdf"               , emit: pca_elbow
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    feature_selection_dimensionality_red.py -ad $input_h5mu -pcs $n_pcs -nn $n_neighbors -md $min_dist

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        feature_selection_dimensionality_red.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix.hvg.h5mu
    touch umap_coordinates.csv
    touch umap_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        feature_selection_dimensionality_red.py --version >> versions.yml
    END_VERSIONS
    """
}
