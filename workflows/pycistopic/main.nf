include { makePseudobulk } from '../../modules/pycistopic/main'
include { peakCalling } from '../../modules/pycistopic/main'
include { inferConsensus } from '../../modules/pycistopic/main'

workflow  PYCISTOPIC {
    take:
        sample_table
        celltype_annotation
        chromsizes
        blacklist
    main:
        // make pseudobulk for each sample
        pseudobulk = makePseudobulk(
            sample_table,
            celltype_annotation,
            chromsizes
        )

        // perform peak calling for pseudobulks
        narrow_peaks = peakCalling(pseudobulk.fragments)

        // get consensus peaks
        consensus = inferConsensus(
            narrow_peaks,
            chromsizes,
            blacklist
        )

        consensus.view()
}