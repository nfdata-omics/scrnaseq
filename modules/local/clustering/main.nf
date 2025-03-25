process CLUSTERING  {
    tag "$meta.id"
    label 'process_single'

    container = 'docker.io/nfdata/sc_rnaseq:v1.0.1'

    input:
    tuple val(meta), path(input_h5mu)

    output:
    //path "combined_matrix_C.h5ad", emit: h5ad
    tuple val(meta), path("*.clustered.h5mu"), emit: h5mu
    //tuple val(run_id), path("combined_matrix_C.h5ad"), emit: h5ad
    //path "ranked_genes.xlsx", emit: DEG
    path "Metadata_final.csv", emit: metadata_final
    path "Leiden_clustering.png", emit: clusters
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    clustering.py -ad $input_h5mu
    
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        clustering.py --version >> versions.yml
    END_VERSIONS 
    """

    stub:
    """
    touch matrix_clustered.h5ad
    touch Metadata_final.csv
    touch Leiden_clustering.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustering.py --version >> versions.yml
    END_VERSIONS
    """
    
    
}