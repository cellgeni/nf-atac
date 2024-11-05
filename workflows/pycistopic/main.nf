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
        cisTopicObjectFlag
    main:
        // Split celltypes and filter sample table from duplicates
        SplitCellTypeAnnotation(sample_table, celltype_annotation)

        // Get splited celltype files
        celltypes = SplitCellTypeAnnotation.out.celltypes

        // Get fragments and barcode metrics paths from filtered sample table
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
        // Conver channgel to list:
        // 1) Collect everything to List with elements [sample_id, fragments_path, fragments_idx_path, barcode_metrics_path]
        // 2) Transpose List -> Channel(sample_id_list, fragments_path_list, fragments_idx_path_list, barcode_metrics_path_list)
        // 3) Convert Channel ->  List[sample_id_list, fragments_path_list, fragments_idx_path_list, barcode_metrics_path_list]
        fragments_list = fragments.collect(flat: false)
                                  .transpose()
                                  .toList()

        // Make pseudobulk for each sample
        pseudobulk = MakePseudobulk(
            fragments_list,
            celltypes,
            chromsizes
        )

        // Perform peak calling for pseudobulks and collect all files
        narrow_peaks = PeakCalling(pseudobulk.fragments)
        narrow_peaks = narrow_peaks.collect()

        // Get consensus peaks
        consensus = InferConsensus(narrow_peaks, chromsizes, blacklist).bed.collect()

        // Perform QC
        fragments_consensus_qc = QualityControl(fragments, consensus, tss_bed)

        // Create cisTopic object
        if (cisTopicObjectFlag) {
            CreateCisTopicObject(fragments_consensus_qc, blacklist)
        }
}