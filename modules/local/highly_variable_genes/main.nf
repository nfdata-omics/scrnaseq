process HIGHLY_VARIABLE_GENES  {
    tag "$meta.id"
    label 'process_single'

    container = 'docker.io/nfdata/sc_rnaseq:v1.0.0'

    input:
    tuple val(meta), path(input_h5ad)

    output:
    tuple val(meta), path("*.hvg.h5ad") , emit: h5ad
    path "UMAP_coordinates.csv", emit: UMAP
    path "UMAP_plot.png", emit: graph_UMAP
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    feature_selection_dimensionality_red.py -ad $input_h5ad

    echo "" >> versions.yml
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
    END_VERSIONS
    feature_selection_dimensionality_red.py --version >> versions.yml

    """

    stub:
    """
    touch matrix.hvg.h5ad
    touch UMAP_coordinates.csv
    touch UMAP_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    feature_selection_dimensionality_red.py --version >> versions.yml

    """
}