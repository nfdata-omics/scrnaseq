process DOUBLETS_ATAC  {
    tag "$meta.id"
    label 'process_single'

    container = 'quay.io/biocontainers/snapatac2:2.8.0--py311h284d45d_1'
    

    input:
    tuple val(meta), path (input_h5ad)
    
    
    output:
    tuple val(meta), path("*.doublets_atac.h5ad"), emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    //doublets_atac.py  -ad $input_h5ad  -id ${meta.collect{ it.id }.join(' ')} 
    
    
    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export XDG_CACHE_HOME=/tmp

    doublets_atac.py  -ad $input_h5ad
    
    
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        doublets_atac.py --version >> versions.yml
    END_VERSIONS
    """
    
    stub:
    """
    touch matrix.doublets_atac.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        doublets_atac.py --version >> versions.yml
    END_VERSIONS
    """
}
