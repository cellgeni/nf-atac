include { SplitCellTypeAnnotation } from '../../modules/pycistopic/main'
include { MakePseudobulk } from '../../modules/pycistopic/main'
include { PeakCalling } from '../../modules/pycistopic/main'
include { CollectPeakMetadata } from '../../modules/pycistopic/main'
include { InferConsensus } from '../../modules/pycistopic/main'
include { QualityControl } from '../../modules/pycistopic/main'
include { CreatePythonObject } from '../../modules/pycistopic/main'
include { CombinePythonObject } from '../../modules/pycistopic/main'
include { CISTOPIC_COUNTFRAGMENTS } from '../../modules/local/cistopic/countfragments'
include { CISTOPIC_SPLITANNOTATION } from '../../modules/local/cistopic/splitannotation'
include { CISTOPIC_PSEUDOBULK } from '../../modules/local/cistopic/pseudobulk'


workflow PEAKCALLING {
    take:
    sample_table
    celltype_annotation
    chromsizes
    fragments_filename
    narrowPeaks_dir
    
    main:
    // STEP1: Get barcode metrics
    // Read sample table and get barcode metrics path
    barcode_metrics = Channel.fromPath(sample_table, checkIfExists: true)
        // Read sample table
        .splitCsv(header: true)
        // Get barcode metrics path
        .map{ row ->
            def barcode_metrics = file( "${row.filedir}/per_barcode_metrics.csv" )
            tuple( [id: row.sample_id], barcode_metrics )
        }
        // Check if barcode metrics file exists
        .branch { meta, barcode_metrics ->
            has_metrics: barcode_metrics.exists()
                return tuple( meta, barcode_metrics)
            no_metrics: !barcode_metrics.exists()
                log.info "No barcode metrics found for sample ${meta.id} at path: ${barcode_metrics}. Calculating..."
                def filedir = barcode_metrics.getParent()
                def fragments_path = file( "${filedir}/${fragments_filename}" )
                // check if fragments file exists
                if ( ! fragments_path.exists() ) {
                    error("No fragments file found for sample ${meta.id} at path: ${fragments_path}. Please check your sample table.")
                }
                return tuple( meta, fragments_path )
        }

    // Count barcode fragments if no barcode metrics found
    CISTOPIC_COUNTFRAGMENTS(barcode_metrics.no_metrics)

    // STEP2: Split celltypes and filter sample table from duplicates
    metrics = barcode_metrics.has_metrics
        .mix(CISTOPIC_COUNTFRAGMENTS.out.csv)
        .collect(flat: false)
        .transpose()
        .toSortedList()
    
    metrics.view()
    CISTOPIC_SPLITANNOTATION(metrics, celltype_annotation, sample_table)

    // Get splited celltype files and combine them with fragment counts
    celltype_fragments = CISTOPIC_SPLITANNOTATION.out.celltype_fragments.splitCsv()
    celltypes = CISTOPIC_SPLITANNOTATION.out.celltypes
                                            .map{ celltype_path -> [ celltype_path.baseName, celltype_path ] }
                                            .transpose()
                                            .combine(celltype_fragments, by: 0)
                                            .map { name, path, fragments -> tuple( [id: name, fragments: fragments], path ) }
    
    // Get fragments files, indexes and fragment counts for each sample (the same as for barcode metrics)
    fragments = CISTOPIC_SPLITANNOTATION.out
                                        .sample_table
                                        .splitCsv(skip: 1)
                                        .map{ sample_id, cellranger_arc_output, _fragment_counts -> 
                                        [
                                            sample_id,
                                            file( "${cellranger_arc_output}/*fragments.tsv.gz" ),
                                            file( "${cellranger_arc_output}/*fragments.tsv.gz.tbi" )
                                        ]
                                        }
                                        .collect(flat: false)
                                        .transpose()
                                        .toList()
    // Make pseudobulk for each sample
    CISTOPIC_PSEUDOBULK(
        fragments,
        celltypes,
        CISTOPIC_SPLITANNOTATION.out.fragments_celltype_x_sample,
        chromsizes
    )

    // Perform peak calling for pseudobulks and collect all files
    PeakCalling(MakePseudobulk.out.pseudobulk_fragments, narrowPeaks_dir)
    narrow_peaks_list = PeakCalling.out.collect(flat: false).transpose().toList()

    CollectPeakMetadata(narrow_peaks_list, narrowPeaks_dir)
    
    emit:
    sample_table = CISTOPIC_SPLITANNOTATION.out.sample_table
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
        CreatePythonObject.out.toList()
                           .branch {
                                it ->
                                combine_objects: it.size() > 1
                                sample: true
                           }
                           .set { objects }

        // Combine cisTopic objects
        CombinePythonObject(objects.combine_objects.transpose().toList())
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
        narrowPeaks_dir
    main:
        // Create pseudobulk, call peaks and update sample table
        if ( callPeaksFlag ) {
            PEAKCALLING(
                sample_table,
                celltypes,
                chromsizes,
                fragments_filename,
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