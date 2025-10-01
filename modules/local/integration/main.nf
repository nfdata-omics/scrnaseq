process INTEGRATION {
    tag "$meta.id"
    label 'process_medium'

    container = 'docker.io/nfdata/muon-sc_rnaseq:v1.0.2'

    input:
    tuple val(meta), path(input_h5mu)

    output:
    tuple val(meta), path("*.integrated.h5mu"), emit: h5mu
    path "versions.yml",  emit: versions

    
    path "Harmony_UMAP_coordinates_GEX.csv", emit: csv_harmony
    path "Harmony-corrected_UMAP_plot_*.png", emit: graph_UMAP_integrated
    

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    # Verify the command that will run
    integration.py -ad $input_h5mu
    

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        integration.py --version >> versions.yml
    END_VERSIONS
    
    
    """
    stub:
    """
    touch matrix.integrated.h5mu
    touch Harmony_UMAP_coordinates_GEX.csv
    touch Harmony-corrected_UMAP_plot_GEX.png


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integration.py --version >> versions.yml
    END_VERSIONS
    
    
    """
}
