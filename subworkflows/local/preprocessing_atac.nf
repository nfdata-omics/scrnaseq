/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { QUALITY_FILTERING_ATAC } from '../../modules/local/qualitymetrics_atac'
include { DIMENSIONALITY_REDUCTION_ATAC } from '../../modules/local/dimensionality_reduction_atac'
include { PEAK_CALLING } from '../../modules/local/peak_calling'

workflow ATAC_PREPROCESSING {

    take:
    fragments
    fragments_index
    tss_threshold
    min_fragments_counts
    max_fragments_counts
    n_features_atac
    frac_dup
    peaks_frac
    n_comps
    n_neighbors
    n_clusters
    blacklist_path
    cell_annotation_meta_ch



    main:
        ch_versions = Channel.empty()


        QUALITY_FILTERING_ATAC (
            fragments,
            fragments_index,
            tss_threshold,
            min_fragments_counts,
            max_fragments_counts,
            blacklist_path

        )
        ch_versions = ch_versions.mix(QUALITY_FILTERING_ATAC.out.versions)


        DIMENSIONALITY_REDUCTION_ATAC (
            QUALITY_FILTERING_ATAC.out.h5ad,
            params.n_features_atac,
            params.frac_dup,
            params.peaks_frac,
            params.n_comps_atac,
            params.n_neighbors_atac,
            params.n_clusters_atac,
            blacklist_path
        )
        ch_versions = ch_versions.mix(DIMENSIONALITY_REDUCTION_ATAC.out.versions)


        PEAK_CALLING(
            DIMENSIONALITY_REDUCTION_ATAC.out.h5ad,
            cell_annotation_meta_ch
        )
        ch_versions = ch_versions.mix(PEAK_CALLING.out.versions)


    emit:
    ch_versions
    h5ad = PEAK_CALLING.out.h5ad_peak



}
