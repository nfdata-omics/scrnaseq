process QUALITY_FILTERING  {
    tag "$meta.id"
    label 'process_single'

    container = 'docker.io/nfdata/sc_atacseq:v1.0.0'
    

    input:
    tuple val(meta), path(input_h5mu)
    tuple val(meta), path (doublets_csv)
    val mt_threshold
    
    output:
    tuple val(meta), path("*.filtered.h5mu"), emit: h5mu
    path "Cells_before_filtering.png", emit: cells_before_filtering, optional: true
    path "Cells_after_filtering.png", emit: cells_after_filtering, optional: true
    path "QC_Density_*.png", emit: qc_density, optional: true
    path "QC_Density_MT-Ribo*.png", emit: qc_density_mito, optional: true
    path "ADTs_Distribution_*.png", emit: adts_distribution, optional: true
    path "Outliers_*.png", emit: outliers, optional: true
    path  "summary_qualitycontrol.csv", emit: summary_qualitycontrol, optional: true
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    
    script:
    
    
    def d = doublets_csv ? "-d $doublets_csv" : ''

    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    

    qualitymetricsfilters.py -ad $input_h5mu $d -mt $mt_threshold
    
    
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        qualitymetricsfilters.py --version >> versions.yml
    END_VERSIONS
    """
    
    stub:
    """
    touch matrix.filtered.h5mu
    touch cells_before_filtering.png
    touch cells_after_filtering.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualitymetricsfilters.py --version >> versions.yml
    END_VERSIONS
    """
}
