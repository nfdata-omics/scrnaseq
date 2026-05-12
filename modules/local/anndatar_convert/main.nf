process ANNDATAR_CONVERT {

    //
    // This module uses the anndata R package to convert h5ad files in different formats
    //

    tag "${meta.id}"

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7f7e0e55a38db6b814fc1d86cca65fbfa60d6ce5e6b2ead0e5db49cd2d83b816/data' :
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-singlecellexperiment_r-seuratobject:603afb8a5c60f65f' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.id}_${meta.input_type}_matrix*.rds"), emit: rds
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'anndatar_convert.R'

    stub:
    """
    touch ${meta.id}_${meta.input_type}_matrix.Rds
    touch versions.yml
    """
}
