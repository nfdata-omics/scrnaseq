process CELL_ANNOTATION_LLM  {
    tag "$meta.id"
    label 'process_high'


    container 'docker.io/nfdata/mllmcelltype:v2.0.4'


    input:
    tuple val(meta), path(input_h5mu)
    val(resolutions)
    val(species)
    val(tissue)

    output:
    tuple val(meta), path("*.annotated.h5mu"), emit: h5mu
    path "umap_*.pdf"                        , emit: graph_umap
    path "annotated_clusters_*.csv"          , emit: annotated_clusters     , optional: true
    path "parameters_*.txt"                  , emit: parameters_txt         , optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export CELLTYPIST_FOLDER=/tmp

    cat > cellannotation_llm.py <<EOF
${file("${moduleDir}/cellannotation_llm.py").text}
EOF

    python3 cellannotation_llm.py \
        --input $input_h5mu \
        --species $species \
        --tissue $tissue \
        --n_genes 10 \
        --pval 0.01 \
        --lfc 1.5 \
        --resolutions $resolutions \
        --umap_embedding X_umap

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
END_VERSIONS
    python3 cellannotation_llm.py --versions-dict "${task.process}" >> versions.yml
    """

    stub:
    """
    touch matrix.annotated.h5mu
    touch umap_leiden_0.4_claude-haiku-4.5_2026.pdf
    touch parameters_2026.txt
    touch annotated_clusters_2026.csv

    cat > cellannotation_llm.py <<EOF
${file("${moduleDir}/cellannotation_llm.py").text}
EOF

    python3 cellannotation_llm.py --versions-dict "${task.process}" > versions.yml
    """
}
