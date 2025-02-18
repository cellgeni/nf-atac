include { SplitCellTypeAnnotation } from '../../modules/pycistopic/main'
include { MakePseudobulk } from '../../modules/pycistopic/main'
include { PeakCalling } from '../../modules/pycistopic/main'
include { CollectPeakMetadata } from '../../modules/pycistopic/main'
include { InferConsensus } from '../../modules/pycistopic/main'
include { QualityControl } from '../../modules/pycistopic/main'
include { CreatePythonObject } from '../../modules/pycistopic/main'
include { CombinePythonObject } from '../../modules/pycistopic/main'


workflow PEAKCALLING {
    take:
        sample_table
        celltype_annotation
        chromsizes
        fragments_filename
        barcode_metrics_filename
        narrowPeaks_dir
    main:
        // Get barcode metrics and convert channel to list:
        // 1) Get barcode metrics paths from sample table
        // 2) Collect everything to List with elements [sample_id, barcode_metrics_path]
        // 3) Transpose List -> Channel(List[sample_id], List[barcode_metrics_path])
        // 4) Convert Channel ->  List(List[sample_id], List[barcode_metrics_path])
        barcode_metrics = Channel.fromPath(sample_table, checkIfExists: true)
                                 .splitCsv(skip: 1)
                                 .map{ sample_id, cellranger_arc_output -> 
                                 [
                                    sample_id,
                                    file( "${cellranger_arc_output}/${barcode_metrics_filename}" ),
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
                                           .map{ sample_id, cellranger_arc_output, _fragment_counts -> 
                                           [
                                                sample_id,
                                                file( "${cellranger_arc_output}/${fragments_filename}" ),
                                                file( "${cellranger_arc_output}/${fragments_filename}.tbi" )
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
        PeakCalling(MakePseudobulk.out.pseudobulk_fragments, narrowPeaks_dir)
        narrow_peaks_list = PeakCalling.out.collect(flat: false).transpose().toList()

        CollectPeakMetadata(narrow_peaks_list, narrowPeaks_dir)
    emit:
        sample_table = SplitCellTypeAnnotation.out.sample_table
        peak_metadata = PeakCalling.out
}


workflow INFERPEAKS {
    take:
        peak_metadata
        sample_table
        chromsizes
        blacklist
        tss_bed
        fragments_filename
    main:
        // Get fragment paths from sample table
        fragments = sample_table.splitCsv(skip: 1)
                                .map{ sample_id, cellranger_arc_output, fragments_num -> 
                                [
                                    sample_id,
                                    file( "${cellranger_arc_output}/${fragments_filename}" ),
                                    file( "${cellranger_arc_output}/${fragments_filename}.tbi" ),
                                    fragments_num
                                ]
                                }

        // Get .narrowPeak paths
        narrow_peak_paths = peak_metadata.map{ _celltype, _fragments, _large_peaks, _all_peaks, path -> path}.collect()

        // Get consensus peaks
        consensus = InferConsensus(narrow_peak_paths, chromsizes, blacklist).bed.collect()

        // Perform QC
        fragments_consensus_qc = QualityControl(fragments, consensus, tss_bed)

        // Create cisTopic object
        CreatePythonObject(fragments_consensus_qc, blacklist)
        objects = CreatePythonObject.out.toList().filter{ it -> it.size() > 1 }.transpose().toList()

        // Combine cisTopic objects
        CombinePythonObject(objects)
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
        fragments_filename
        barcode_metrics_filename
        narrowPeaks_dir
    main:
        // Create pseudobulk, call peaks and update sample table
        if ( callPeaksFlag ) {
            PEAKCALLING(
                sample_table,
                celltypes,
                chromsizes,
                fragments_filename,
                barcode_metrics_filename,
                narrowPeaks_dir
            )
            // get pseudobulk peaks and updated sample table
            peak_metadata = PEAKCALLING.out.peak_metadata
            sample_table = PEAKCALLING.out.sample_table
        } else {
            sample_table = Channel.fromPath(sample_table, checkIfExists: true)
            peak_metadata = Channel.fromPath(celltypes, checkIfExists: true).splitCsv(skip:1, sep:'\t')
        }

        if ( inferConsensusFlag ) {
            INFERPEAKS(
                peak_metadata,
                sample_table,
                chromsizes,
                blacklist,
                tss_bed,
                fragments_filename
            )
        }
}