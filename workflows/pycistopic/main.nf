include { MakePseudobulk } from '../../modules/pycistopic/main'
include { PeakCalling } from '../../modules/pycistopic/main'
include { InferConsensus } from '../../modules/pycistopic/main'
include { QualityControl } from '../../modules/pycistopic/main'
include { CreateCisTopicObject } from '../../modules/pycistopic/main'

workflow  PYCISTOPIC {
    take:
        sample_table
        celltype_annotation
        chromsizes
        blacklist
        tss_bed
    main:
        // Make pseudobulk for each sample
        pseudobulk = MakePseudobulk(
            sample_table,
            celltype_annotation,
            chromsizes
        )

        // Perform peak calling for pseudobulks
        narrow_peaks = PeakCalling(pseudobulk.fragments)

        // Get consensus peaks
        consensus = InferConsensus(
            narrow_peaks,
            chromsizes,
            blacklist
        )

        // Perform QC
        fragments_consensus = sample_table.join(consensus, failOnDuplicate: true)
        fragments_consensus_qc = QualityControl(fragments_consensus, tss_bed)

        // Create cisTopic object
        CreateCisTopicObject(fragments_consensus_qc, blacklist)
}