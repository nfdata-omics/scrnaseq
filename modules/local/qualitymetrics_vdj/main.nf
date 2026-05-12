process QUALITY_METRICS_VDJ   {
    tag "$meta.id"
    label 'process_medium'

    container = 'quay.io/biocontainers/scirpy:0.20.1--pyhdfd78af_0'


    input:
    tuple val(meta), path(input_h5mu)


    output:
    tuple val(meta), path("*.qc_vdj.h5mu"),     emit: h5mu
    path "VDJ_abundance_all_samples.csv",       emit: vdj_abundance,       optional: true
    path "vdj_receptor_type.pdf",               emit: receptor_type,       optional: true
    path "vdj_receptor_subtype.pdf",            emit: receptor_subtype,    optional: true
    path "vdj_chain_pairing.pdf",               emit: chain,               optional: true
    path "chain_pairing_stats_all_samples.csv", emit: chain_pairing_stats, optional: true
    path "vdj_gex_match_all_samples.csv",       emit: gex_match,           optional: true
    path "ChainPairing_donor-*.csv",            emit: chain_pairing_donor, optional: true
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:


    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp



    qualitymetrics_vdj.py -ad $input_h5mu


    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        qualitymetrics_vdj.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix.qc_vdj.h5mu


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualitymetrics_vdj.py --version >> versions.yml
    END_VERSIONS
    """
}
