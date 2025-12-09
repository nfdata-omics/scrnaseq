process INTEGRATION {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/nfdata/muon-sc_rnaseq:v1.0.2'

    input:
    tuple val(meta), path(input_h5mu)
    val n_neighbors_harmony
    val min_dist_harmony

    output:
    tuple val(meta), path("*.integrated.h5mu")           , emit: h5mu
    path "versions.yml"                                  , emit: versions
    path "Harmony_UMAP_coordinates_GEX.csv"              , emit: csv_harmony
    path "Harmony_corrected_UMAP_plot_GEX.pdf"           , emit: graph_UMAP_integrated
    path "Harmony_corrected_UMAP_plot_GEX_QC.pdf"        , emit: graph_UMAP_integrated_QC
    path "Harmony_corrected_UMAP_plot_GEX_phase.pdf"     , emit: graph_UMAP_integrated_phase, optional: true
    path "Harmony_corrected_UMAP_plot_GEX_celltypist.pdf", emit: graph_UMAP_integrated_celltypist, optional: true
    path "Harmony_corrected_UMAP_plot_GEX_metadata.pdf"  , emit: graph_UMAP_integrated_metadata, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp


    integration.py -ad $input_h5mu -nnh $n_neighbors_harmony -mdh $min_dist_harmony


    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        integration.py --version >> versions.yml
    END_VERSIONS


    """
    stub:
    """
    touch matrix.integrated.h5mu
    touch Harmony_UMAP_coordinates_GEX.csv
    touch Harmony-corrected_UMAP_plot_GEX.pdf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integration.py --version >> versions.yml
    END_VERSIONS


    """
}
