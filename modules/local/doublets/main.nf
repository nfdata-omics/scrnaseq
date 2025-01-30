process DOUBLETS  {
    tag "$meta.id"
    label 'process_single'
    
    container = 'docker.io/nfdata/sc-rnaseq-r:v1.0.0'

    input:
    tuple val(meta), path(input_sce)
    

    output:
    path "doublets_score.csv", emit: doublets
    

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    doublets.R ${input_sce[0]}
    
    
    """
    
    stub:

    """
    touch doublets_score.csv
    
    
    
    """
}