include { makePseudobulk } from '../../modules/pycistopic/main'

workflow  PYCISTOPIC {
    take:
        sample_table
        celltype_annotation
        chromsizes
    main:
        // make pseudobulk for each sample
        pseudobulk = makePseudobulk(
            sample_table,
            celltype_annotation,
            chromsizes
        )
        pseudobulk.bed.view()
        pseudobulk.bigwig.view()
}