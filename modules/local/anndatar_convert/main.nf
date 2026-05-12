process ANNDATAR_CONVERT {

    //
    // This module uses the anndata R package to convert h5ad files in different formats
    //

    tag "${meta.id}"

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b6f6a502639a53b3b7b078b1b250c0c4a953cfd5508feac34c4d8185e3b2de24/data' :
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-singlecellexperiment_r-hdf5r_r-seurat:67d97559705c8ef0' }"

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
