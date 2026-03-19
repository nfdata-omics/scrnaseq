include { PSEUDOBULK            } from '../../modules/local/pseudobulk/main.nf'

workflow PSEUDOBULK_ANALYSIS {

    take:

    h5mu
    resolution
    group_column
    comparisons
    
    main:

    h5mu
        .combine(resolution)
        .combine(group_column)
        .combine(comparisons)
         .map { meta, h5mu, resolution, group_column, comparisons -> tuple(meta, h5mu, group_column, resolution, comparisons) }
         .set { pseudobulk_inputs }

    PSEUDOBULK(
        pseudobulk_inputs
    )

    emit:

    pseudobulk_results = PSEUDOBULK.out.pseudobulk_results
    versions = PSEUDOBULK.out.versions

}