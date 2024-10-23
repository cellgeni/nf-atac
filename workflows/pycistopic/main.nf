include { makePseudobulk } from '../../modules/pycistopic/main'
include { peakCalling } from '../../modules/pycistopic/main'

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

        peakCalling = peakCalling(pseudobulk.bed)
        peakCalling.view()
}