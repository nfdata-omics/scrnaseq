process DIMENSIONALITY_REDUCTION_ATAC  {
    tag "$meta.id"
    label 'process_medium'

    container = 'quay.io/biocontainers/snapatac2:2.8.0--py311h284d45d_1'
    

    input:
    tuple val(meta), path (input_h5ad)
    
    
    output:
    tuple val(meta), path("*.dimred_atac.h5ad"), emit: h5ad
    path "umap_ATAC.png", emit: umap_atac, optional: true
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    
    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export XDG_CACHE_HOME=/tmp

    dimensionalityreduction_atac.py  -ad $input_h5ad
    
    
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        dimensionalityreduction_atac.py --version >> versions.yml
    END_VERSIONS
    """
    
    stub:
    """
    touch matrix.dimred_atac.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dimensionalityreduction_atac.py --version >> versions.yml
    END_VERSIONS
    """
}
