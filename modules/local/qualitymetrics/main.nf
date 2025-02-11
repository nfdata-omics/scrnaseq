process QUALITY_FILTERING  {
    tag "$meta.id"
    label 'process_single'

    container = 'docker.io/nfdata/sc_rnaseq:v1.0.0'

    input:
    tuple val(meta), path(input_h5mu)
    path doublets_csv
    val mt_threshold
    
    output:
    tuple val(meta), path("*.filtered.h5ad"), emit: h5mu
    path "Cells_before_filtering.png", emit: cells_before_filtering
    path "Cells_after_filtering.png", emit: cells_after_filtering
    path "QC_Density_*.png", emit: qc_density
    path "QC_Density_MT-Ribo*.png", emit: qc_density_mito
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    qualitymetricsfilters.py -ad $input_h5mu -d $doublets_csv -f $mt_threshold
    
    echo "" >> versions.yml
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
    END_VERSIONS
    qualitymetricsfilters.py --version >> versions.yml
    
    """
    
    stub:
    """
    touch matrix.filtered.h5ad
    touch cells_before_filtering.png
    touch cells_after_filtering.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    qualitymetricsfilters.py --version >> versions.yml
    
    """
}