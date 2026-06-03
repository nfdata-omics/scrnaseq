process DESEQ2_COMPARE {
    tag "${meta.id}_${comparison}"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(model), val(comparison)
    val fdr_threshold

    output:
    tuple val(meta), path("deseq2_toptable.*.txt"), emit: dge
    tuple val(meta), path("deseq2_summary.*.txt") , emit: summary
    tuple val(meta), path("pvalue_hist.*.pdf")    , emit: hist_pval
    tuple val(meta), path("volcano.*.pdf")        , emit: volcano
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DESEQ2_COMPARE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    cat > dge_deseq2_results.R <<'RSCRIPT'
${file("${moduleDir}/dge_deseq2_results.R").text}
RSCRIPT

    Rscript dge_deseq2_results.R $args --FDR $fdr_threshold $model "$comparison"
    cat <(echo -e '#DGE\t${meta.id}\t${comparison}') deseq2_toptable.*.txt > out_tmp
    j=\$(basename deseq2_toptable.*.txt)
    mv out_tmp "\$j"

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
END_VERSIONS
    Rscript dge_deseq2_results.R --version >> versions.yml
    """

    stub:
    """
    touch deseq2_toptable.${meta.id}_${comparison}.txt
    touch deseq2_summary.${meta.id}_${comparison}.txt
    touch pvalue_hist.${meta.id}_${comparison}.pdf
    touch volcano.${meta.id}_${comparison}.pdf

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
END_VERSIONS
    Rscript dge_deseq2_results.R --version >> versions.yml
    """
}
