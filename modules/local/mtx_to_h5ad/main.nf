process MTX_TO_H5AD {

    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::scanpy==1.10.2 conda-forge::python-igraph conda-forge::leidenalg"
    container "community.wave.seqera.io/library/scanpy:1.10.2--e83da2205b92a538"

    input:
    tuple val(meta), path(matrix_paths)

    output:
    tuple val(meta), path("${meta.id}_${meta.input_type}_matrix.h5ad"), emit: h5ad
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def matrix_input = selectMatrixInput(matrix_paths, meta)

    template "mtx_to_h5ad_cellranger.py"

    stub:
    """
    touch ${meta.id}_${meta.input_type}_matrix.h5ad
    touch versions.yml
    """
}

/*
 * Cell Ranger modules currently expose all files below their output directory,
 * whereas samplesheet restart entries contain one logical matrix artifact.
 * Prefer the H5 representation when both H5 and MEX are present in the legacy
 * aligner output. A supplied MEX directory is used directly.
 */
def selectMatrixInput(matrix_paths, meta) {
    def paths = matrix_paths instanceof Collection ? matrix_paths as List : [matrix_paths]
    def h5_paths = paths.findAll { path ->
        path.name.toLowerCase().endsWith('.h5')
    }

    if (h5_paths.size() == 1) {
        return h5_paths.first()
    }
    if (h5_paths.size() > 1) {
        error("MTX_TO_H5AD received multiple H5 matrices for sample '${meta.id}' (${meta.input_type}): ${h5_paths*.name.join(', ')}")
    }

    def mex_paths = paths.findAll { path ->
        isMexDirectory(path)
    }
    if (mex_paths.size() == 1) {
        return mex_paths.first()
    }
    if (mex_paths.size() > 1) {
        error("MTX_TO_H5AD received multiple MEX directories for sample '${meta.id}' (${meta.input_type}): ${mex_paths*.name.join(', ')}")
    }

    error("MTX_TO_H5AD requires one H5 matrix or one complete MEX directory for sample '${meta.id}' (${meta.input_type}).")
}

def isMexDirectory(path) {
    def input_path = path as java.nio.file.Path
    if (!java.nio.file.Files.isDirectory(input_path)) {
        return false
    }

    def required_components = [
        ['matrix.mtx', 'matrix.mtx.gz'],
        ['barcodes.tsv', 'barcodes.tsv.gz'],
        ['features.tsv', 'features.tsv.gz']
    ]
    return required_components.every { alternatives ->
        alternatives.any { filename ->
            java.nio.file.Files.isRegularFile(input_path.resolve(filename))
        }
    }
}
