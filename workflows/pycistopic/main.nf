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
        // Split celltypes and filter sample table from duplicates
        SplitCellTypeAnnotation(sample_table, celltype_annotation)

        // Get splited celltype files
        celltypes = SplitCellTypeAnnotation.out.celltypes

        // Get fragments and barcode metrics from filtered sample table. Steps' description:
        // 1) split sample_table in rows (skipping header)
        // 2) convert cellranger-arc output dir to fragments paths
        // 3) convert everything to List with elements [sample_id, fragments_path, fragments_idx_path, barcode_metrics_path]
        // 4) Transpose List -> Channel(sample_id_list, fragments_path_list, fragments_idx_path_list, barcode_metrics_path_list)
        // 5) Convert Channel ->  List[sample_id_list, fragments_path_list, fragments_idx_path_list, barcode_metrics_path_list]
        fragments = SplitCellTypeAnnotation.out
                                          .sample_table
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

        //Make pseudobulk for each sample
        pseudobulk = MakePseudobulk(
            fragments,
            celltypes,
            chromsizes
        )

        pseudobulk.fragments.view()

        // Perform peak calling for pseudobulks
        narrow_peaks = PeakCalling(pseudobulk.fragments)
        narrow_peaks = narrow_peaks.collect()

        // // Get consensus peaks
        consensus = InferConsensus(
            narrow_peaks,
            chromsizes,
            blacklist
        )

        // // Perform QC
        // fragments_consensus = sample_table.join(consensus, failOnDuplicate: true)
        // fragments_consensus_qc = QualityControl(fragments_consensus, tss_bed)

        // // Create cisTopic object
        // CreateCisTopicObject(fragments_consensus_qc, blacklist)
}