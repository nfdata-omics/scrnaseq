process PSEUDOBULK {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/nfdata/muon-sc_rnaseq:v1.0.4'

    input:
    tuple val(meta), path(h5mu), val(group_column), val(resolution), val(comparison)

    output:
    path "pseudobulk/*", emit: pseudobulk_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    cat > pseudobulk.py <<EOF
${file("${moduleDir}/pseudobulk.py").text}
EOF

    python3 pseudobulk.py \
        --mdata $h5mu \
        --resolution $resolution \
        --group_column $group_column \
        --comparisons $comparison

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
END_VERSIONS
    python3 pseudobulk.py --versions-dict "${task.process}" >> versions.yml
    """

    stub:
    """
    mkdir pseudobulk/

    cat > pseudobulk.py <<EOF
${file("${moduleDir}/pseudobulk.py").text}
EOF

    python3 pseudobulk.py --versions-dict "${task.process}" > versions.yml
    """
}