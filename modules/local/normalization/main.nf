process NORMALIZATION   {
    tag "$meta.id"
    label 'process_medium'

    container = 'docker.io/nfdata/muon-sc_rnaseq:v1.0.1'

    input:
    tuple val(meta), path(input_h5mu)
    tuple val(meta), path(input_raw_h5ad)
    path input_cellcycle_file

    output:
    tuple val(meta), path("*.norm.h5mu"), emit: h5mu
    path "pca_cellcycle_GEX_phase.png", emit: pca_cellcycle_GEX_phase,optional: true
    path "pca_cellcycle_GEX_sample.png", emit: pca_cellcycle_GEX_sample,optional: true
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def rw = input_raw_h5ad ? "-rw $input_raw_h5ad" : ''

    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    normalization.py -ad $input_h5mu $rw --input_cellcycle_file $input_cellcycle_file

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        normalization.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix.norm.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        normalization.py --version >> versions.yml
    END_VERSIONS
    """

}
