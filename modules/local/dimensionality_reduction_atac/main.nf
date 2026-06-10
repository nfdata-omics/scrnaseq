process DIMENSIONALITY_REDUCTION_ATAC  {
    tag "$meta.id"
    label 'process_high_memory'

    container 'docker.io/nfdata/snapatac:v1.0.2'

    input:
    tuple val(meta), path (input_h5ad)
    val  n_features_atac
    val frac_dup
    val peaks_frac
    val n_comps_atac
    val n_neighbors_atac
    val n_clusters_atac
    path blacklist_path
    path cell_counts

    output:
    tuple val(meta), path("*.dimred_atac.h5ad"), emit: h5ad
    path "umap_ATAC_sample_before_Harmony.pdf", emit: umap_sample_before_harmony, optional: true
    path "umap_ATAC_sample_Harmony.pdf", emit: umap_sample_harmony, optional: true
    path "umap_ATAC_all_res_Harmony.pdf", emit: umap_all_res_harmony, optional: true
    path "umap_ATAC_res_*_Harmony.png", emit: umap_individual_harmony, optional: true
    path "umap_ATAC_meta_Harmony.pdf", emit: umap_metadata_harmony, optional: true
    path "Cells_after_filtering_atac.pdf", emit: qc_plots, optional: true
    path "cell_counts.csv", emit: cell_counts, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp
    export XDG_CACHE_HOME=/tmp

    dimensionalityreduction_atac.py  -ad $input_h5ad -f $n_features_atac -fd $frac_dup -pf $peaks_frac -nc $n_comps_atac -nn $n_neighbors_atac -ncl $n_clusters_atac -b $blacklist_path -cc $cell_counts


    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        dimensionalityreduction_atac.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch matrix.dimred_atac.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dimensionalityreduction_atac.py --version >> versions.yml
    END_VERSIONS
    """
}
