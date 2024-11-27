include { SplitCellTypeAnnotation } from '../../modules/pycistopic/main'
include { MakePseudobulk } from '../../modules/pycistopic/main'
include { PeakCalling } from '../../modules/pycistopic/main'
include { InferConsensus } from '../../modules/pycistopic/main'
include { QualityControl } from '../../modules/pycistopic/main'
include { CreateCisTopicObject } from '../../modules/pycistopic/main'


workflow PEAKCALLING {
    take:
        sample_table
        celltype_annotation
        chromsizes
    main:
        // Get barcode metrics and convert channel to list:
        // 1) Get barcode metrics paths from sample table
        // 2) Collect everything to List with elements [sample_id, barcode_metrics_path]
        // 3) Transpose List -> Channel(sample_id_list, barcode_metrics_path_list)
        // 4) Convert Channel ->  List[sample_id_list, barcode_metrics_path_list]
        barcode_metrics = Channel.fromPath(sample_table, checkIfExists: true)
                                 .splitCsv(skip: 1)
                                 .map{ sample_id, cellranger_arc_output -> 
                                 [
                                    sample_id,
                                    file( "${cellranger_arc_output}/${params.barcode_metrics_filename}" ),
                                 ]}
                                 .collect(flat: false)
                                 .transpose()
                                 .toList()
  

        // Split celltypes and filter sample table from duplicates
        SplitCellTypeAnnotation(barcode_metrics, celltype_annotation, sample_table)

        // Get splited celltype files and combine them with fragment counts
        celltype_fragments = SplitCellTypeAnnotation.out.celltype_fragments.splitCsv()
        celltypes = SplitCellTypeAnnotation.out.celltypes
                                               .map{ celltype_path -> [ celltype_path.baseName, celltype_path ] }
                                               .transpose()
                                               .combine(celltype_fragments, by: 0)
        
        // Get fragments files, indexes and fragment counts for each sample (the same as for barcode metrics)
        fragments = SplitCellTypeAnnotation.out
                                           .sample_table
                                           .splitCsv(skip: 1)
                                           .map{ sample_id, cellranger_arc_output, fragment_counts -> 
                                           [
                                                sample_id,
                                                file( "${cellranger_arc_output}/${params.fragments_filename}" ),
                                                file( "${cellranger_arc_output}/${params.fragments_filename}.tbi" )
                                           ]
                                           }
                                           .collect(flat: false)
                                           .transpose()
                                           .toList()
        // Make pseudobulk for each sample
        MakePseudobulk(
            fragments,
            celltypes,
            SplitCellTypeAnnotation.out.fragments_celltype_x_sample,
            chromsizes
        )

        // Perform peak calling for pseudobulks and collect all files
        narrow_peaks = PeakCalling(MakePseudobulk.out.pseudobulk_fragments)
    publish:
        narrow_peaks >> 'narrowPeaks'
    emit:
        narrow_peaks = narrow_peaks
        sample_table = SplitCellTypeAnnotation.out.sample_table
}


workflow INFERPEAKS {
    take:
        narrow_peaks
        sample_table
        chromsizes
        blacklist
        tss_bed
    main:
        // Get fragment paths from sample table
        fragments = Channel.fromPath(sample_table, checkIfExists: true)
                            .splitCsv(skip: 1)
                            .map{ sample_id, cellranger_arc_output, fragments_num -> 
                            [
                                sample_id,
                                file( "${cellranger_arc_output}/${params.fragments_filename}" ),
                                file( "${cellranger_arc_output}/${params.fragments_filename}.tbi" ),
                                fragments_num
                            ]
                            }

        // Get .narrowPeak paths
        narrow_peak_paths = narrow_peaks.map{ celltype, fragments, large_peaks, all_peaks, path -> path}

        // Get consensus peaks
        consensus = InferConsensus(narrow_peak_paths, chromsizes, blacklist).bed.collect()

        // Perform QC
        fragments_consensus_qc = QualityControl(fragments, consensus, tss_bed)

        // Create cisTopic object
        CreateCisTopicObject(fragments_consensus_qc, blacklist)
}


workflow  PYCISTOPIC {
    take:
        sample_table
        celltypes
        chromsizes
        blacklist
        tss_bed
        callPeaksFlag
        inferConsensusFlag
    main:
        // Create pseudobulk, call peaks and update sample table
        if ( callPeaksFlag ) {
            PEAKCALLING(
                sample_table,
                celltypes,
                chromsizes
            )
            // get pseudobulk peaks and updated sample table
            pseudobulk_peaks = PEAKCALLING.out.narrow_peaks
            sample_table = PEAKCALLING.out.sample_table
        } else {
            narrow_peaks = Channel.fromPath(celltypes, checkIfExists: true).splitCsv()
        }

        if ( inferConsensusFlag ) {
            INFERPEAKS(
                narrow_peaks,
                sample_table,
                chromsizes,
                blacklist,
                tss_bed
            )
        }
}