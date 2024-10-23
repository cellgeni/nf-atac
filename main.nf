include { PYCISTOPIC } from './workflows/pycistopic/main'


workflow {
    // Read sample table and convert string with path to Path object
    sample_table = Channel
                        .fromPath(params.sample_table, checkIfExists: true)
                        .splitCsv(skip: 1)
                        .map{ sample_id, fragments_path -> [sample_id, file( fragments_path ), file( "${fragments_path}.tbi" )]}

    // Load celltype annotation and chromsizes
    celltype_annotation = file( params.celltype_annotation )
    chromsizes = file( params.chromsizes )

    // Run PyCistopic pipeline
    PYCISTOPIC(
        sample_table,
        celltype_annotation,
        chromsizes
    )

    // PYCISTOPIC.out.bed.view()
    // PYCISTOPIC.out.bigwig.view()
}