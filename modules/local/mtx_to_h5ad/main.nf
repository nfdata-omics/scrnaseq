process MTX_TO_H5AD {

    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::scanpy==1.10.2 conda-forge::python-igraph conda-forge::leidenalg"
    container "community.wave.seqera.io/library/scanpy:1.10.2--e83da2205b92a538"

    input:
    tuple val(meta), path(matrix_input)

    output:
    tuple val(meta), path("${meta.id}_${meta.input_type}_matrix.h5ad"), emit: h5ad
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "mtx_to_h5ad_cellranger.py"

    stub:
    """
    touch ${meta.id}_${meta.input_type}_matrix.h5ad
    touch versions.yml
    """
}
