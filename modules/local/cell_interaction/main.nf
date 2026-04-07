process CELL_INTERACTION {
    tag "${meta.id}_${resolution}"
    label 'process_single'

    container 'docker.io/nfdata/cell-interaction:v1.7.1'

    input:
    tuple val(meta), path(h5mu), val(method), val(resource), val(resolution)

    output:
    path "*.pdf", emit: plots
    path "*.xlsx", emit: excel
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    cat > cell_interaction.py <<EOF
${file("${moduleDir}/cell_interaction.py").text}
EOF

    python3 cell_interaction.py \
        -ad $h5mu \
        -m $method \
        -resource $resource \
        --results cell_interaction \
        --resolution $resolution

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        python3 cell_interaction.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    mkdir cell_interaction/

    cat > cell_interaction.py <<EOF
${file("${moduleDir}/cell_interaction.py").text}
EOF

    python3 cell_interaction.py --versions-dict "${task.process}" > versions.yml
    """

}
