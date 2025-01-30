process QUALITY_FILTERING  {
    tag "$meta.id"
    label 'process_single'

    container = 'docker.io/nfdata/sc_rnaseq:v1.0.0'

    input:
    tuple val(meta), path(input_h5ad)
    path doublets_csv
    val MT
    
    

    output:
    tuple val(meta), path("*.filtered.h5ad"), emit: h5ad
    path "Cells_before_filtering.png", emit: cells_beforefiltering
    path "Cells_after_filtering.png", emit: cells_afterfiltering
    path "QC_Density_10k_PBMC_5pv2_nextgem_Chromium_X_gex_1_filtered.png", emit: QC_X
    path "QC_Density_10k_PBMC_5pv2_nextgem_Chromium_X_gex_2_filtered.png", emit: QC_Y
    path "QC_Density_MT-Ribo_10k_PBMC_5pv2_nextgem_Chromium_X_gex_1_filtered.png", emit: mitochondrial_X
    path "QC_Density_MT-Ribo_10k_PBMC_5pv2_nextgem_Chromium_X_gex_2_filtered.png", emit: mitochondrial_Y
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    quality_metrics-filters.py -ad $input_h5ad -d $doublets_csv -f $MT
    


    echo "" >> versions.yml
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
    END_VERSIONS
    quality_metrics-filters.py --version >> versions.yml
    
    """
    
    stub:
    """
    touch matrix.filtered.h5ad
    touch Cells_before_filtering.png
    touch Cells_after_filtering.png
    //touch QC_Density_${meta}.png
    //touch QC_Density_MT-Ribo_${meta.original_id}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    quality_metrics-filters.py --version >> versions.yml
    

    """
}