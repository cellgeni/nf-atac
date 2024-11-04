include { SplitCellTypeAnnotation } from '../../modules/pycistopic/main'
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
        fragments = Channel
                        .fromPath(sample_table, checkIfExists: true)
                        .splitCsv(skip: 1)
                        .map{ sample_id, cellranger_arc_output -> 
                                [
                                    sample_id,
                                    file( "${cellranger_arc_output}/${params.fragments_filename}" ),
                                    file( "${cellranger_arc_output}/${params.fragments_filename}.tbi" ),
                                    file( "${cellranger_arc_output}/${params.barcode_metrics_filename}" ),
                                ]
                        }
                        .collect(flat: false)
                        .transpose()
                        .toList()
        // Split celltypes
        celltypes = SplitCellTypeAnnotation(sample_table, celltype_annotation)

        // Make pseudobulk for each sample
        pseudobulk = MakePseudobulk(
            fragments,
            celltypes,
            chromsizes
        )

        pseudobulk.fragments.view()

        // Perform peak calling for pseudobulks
        // narrow_peaks = PeakCalling(pseudobulk.fragments)

        // // Get consensus peaks
        // consensus = InferConsensus(
        //     narrow_peaks,
        //     chromsizes,
        //     blacklist
        // )

        // // Perform QC
        // fragments_consensus = sample_table.join(consensus, failOnDuplicate: true)
        // fragments_consensus_qc = QualityControl(fragments_consensus, tss_bed)

        // // Create cisTopic object
        // CreateCisTopicObject(fragments_consensus_qc, blacklist)
}