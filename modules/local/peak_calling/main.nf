process PEAK_CALLING  {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/nfdata/snapatac:v1.0.0'


    input:
    tuple val(meta), path (input_h5ad)
    path input_meta_file


    output:
    tuple val(meta), path("matrix.tile_atac.h5ad"), emit: h5ad_tile, optional: true
    tuple val(meta), path("matrix.peak_atac.h5ad"), emit: h5ad_peak
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export XDG_CACHE_HOME=/tmp

    peak_calling.py  -ad $input_h5ad -meta $input_meta_file

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        peak_calling.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix.tile_atac.h5ad
    touch matrix.peak_atac.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peak_calling.py --version >> versions.yml
    END_VERSIONS
    """
}
