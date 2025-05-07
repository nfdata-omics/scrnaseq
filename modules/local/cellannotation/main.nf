process CELL_ANNOTATION  {
    tag "$meta.id"
    label 'process_single'


    container = 'docker.io/nfdata/muon-sc_rnaseq:v1.0.3'
    //container = 'quay.io/teichlab/celltypist:latest'

    input:
    tuple val(meta), path(input_h5mu)
    path input_model

    output:
    tuple val(meta), path("*.annotated.h5mu") , emit: h5mu
    path "umap_plot_*.png", emit: graph_umap
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export CELLTYPIST_FOLDER=/tmp

    cellannotation.py -ad $input_h5mu --model-list $input_model

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        cellannotation.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix.annotated.h5mu
    touch umap_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellannotation.py --version >> versions.yml
    END_VERSIONS
    """
}
