include { CISTOPIC_COUNTFRAGMENTS } from '../../modules/local/cistopic/countfragments'
include { CISTOPIC_SPLITANNOTATION } from '../../modules/local/cistopic/splitannotation'
include { CISTOPIC_PSEUDOBULK } from '../../modules/local/cistopic/pseudobulk'
include { CISTOPIC_CALLPEAKS } from '../../modules/local/cistopic/callpeaks'
include { CISTOPIC_INFERCONSENSUS } from '../../modules/local/cistopic/inferconsensus'
include { CISTOPIC_QUALITYCONTROL } from '../../modules/local/cistopic/qualitycontrol'
include { CISTOPIC_CREATEOBJECT } from '../../modules/local/cistopic/createobject'
include { CISTOPIC_COMBINEOBJECTS } from '../../modules/local/cistopic/combineobjects'
include { ANNDATA_CONCAT } from '../../modules/local/anndata/concat'

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
        .collect(flat: false, sort: true)
        .transpose()
        .toList()
    
    CISTOPIC_SPLITANNOTATION(metrics, celltype_annotation, sample_table)

    // STEP3: Make pseudobulks for each celltype
    // Get splited celltype files and combine them with fragment counts
    celltype_fragments = CISTOPIC_SPLITANNOTATION.out.celltype_fragments.splitCsv()
    celltypes = CISTOPIC_SPLITANNOTATION.out.celltypes
        .map{ celltype_path -> [ celltype_path.baseName, celltype_path ] }
        .transpose()
        .combine(celltype_fragments, by: 0)
        .map { name, path, fragments -> tuple( [id: name, fragments: fragments], path ) }

    // Get fragments files, indexes and fragment counts for each sample (the same as for barcode metrics)
    fragments = CISTOPIC_SPLITANNOTATION.out.sample_table
        .splitCsv(header: true)
        .map{ meta -> 
            [
                [id: meta.sample_id, fragments: meta.fragments],
                file( "${meta.filedir}/*fragments.tsv.gz" )[0],
                file( "${meta.filedir}/*fragments.tsv.gz.tbi" )[0]
            ]
        }
        .collect(flat: false, sort: true)
        .transpose()
        .toList()

    // Make pseudobulk for each sample
    CISTOPIC_PSEUDOBULK(
        fragments,
        celltypes,
        chromsizes,
        CISTOPIC_SPLITANNOTATION.out.fragments_celltype_x_sample
    )

    // STEP 4: Perform peak calling for pseudobulks and collect all files
    CISTOPIC_CALLPEAKS(CISTOPIC_PSEUDOBULK.out.tsv)

    // Collect peak metadata
    peaks = CISTOPIC_CALLPEAKS.out.narrowPeak.map { meta, path, large_peaks, all_peaks -> 
            def updated_meta = meta + [ large_peaks: large_peaks, all_peaks: all_peaks ]
            tuple( updated_meta, path )
        }
    output_dir = params.output_dir ? params.output_dir : "."
    peaks.collectFile(
            name: "pseudobulk_peaks.csv",
            storeDir: output_dir,
            newLine: true,
            seed: "celltype,fragments,large_peaks,all_peaks,path",
            sort: { line -> line.split(',')[0] }
        ) { meta, path ->
            "${meta.id},${meta.fragments},${meta.large_peaks},${meta.all_peaks},${path.toString()}"
        }
        .subscribe { __ -> 
                log.info("Pseudobulk peak calling stats are saved to ${output_dir}/pseudobulk_peaks.tsv")
            }
    // Collect versions
    versions = CISTOPIC_COUNTFRAGMENTS.out.versions.first()
        .mix(
            CISTOPIC_SPLITANNOTATION.out.versions,
            CISTOPIC_PSEUDOBULK.out.versions.first(),
            CISTOPIC_CALLPEAKS.out.versions.first()
        )

    emit:
    sample_table = CISTOPIC_SPLITANNOTATION.out.sample_table
    pseudobulk   = CISTOPIC_PSEUDOBULK.out.tsv
    peaks        = peaks
    versions     = versions
}


workflow INFERPEAKS {
    take:
        peaks
        sample_table
        chromsizes
        blacklist
        tss_bed
        fragments_filename
    main:
        // STEP 0: Prepare inputs
        // Get fragment paths from sample table
        fragments = sample_table
            .splitCsv(header: true)
            .map{ meta -> 
                [
                    [id: meta.sample_id, fragments: meta.fragments],
                    file( "${meta.filedir}/*fragments.tsv.gz" )[0],
                    file( "${meta.filedir}/*fragments.tsv.gz.tbi" )[0]
                ]
            }

        // Get peak paths from peaks channel
        narrowPeaks = peaks.toSortedList()
            .map{ list -> 
                def celltype_names = list.collect{ it[0].id }
                def narrowpeak_files = list.collect{ it[1] }
                return [ celltype_names, narrowpeak_files ]
            }
        
        // STEP 1: Get consensus peaks
        CISTOPIC_INFERCONSENSUS(
            narrowPeaks,
            chromsizes,
            blacklist
        )

        consensus = CISTOPIC_INFERCONSENSUS.out.bed.collect()

        // STEP 2: Perform QC
        CISTOPIC_QUALITYCONTROL(fragments, consensus, tss_bed)

        // SPET 3: Create cisTopic object
        CISTOPIC_CREATEOBJECT(
            CISTOPIC_QUALITYCONTROL.out.qc,
            consensus,
            blacklist
        )
        
        cistopic_objects = CISTOPIC_CREATEOBJECT.out.pkl
            .toSortedList()
            .branch {
                it ->
                combine_objects: it.size() > 1
                sample: true
            }
        
        anndata_objects = CISTOPIC_CREATEOBJECT.out.h5ad
            .toSortedList()
            .branch {
                it ->
                combine_objects: it.size() > 1
                sample: true
            }

        // STEP 4: Combine cisTopic objects if needed
        CISTOPIC_COMBINEOBJECTS(cistopic_objects.combine_objects.transpose().toList())
        ANNDATA_CONCAT(anndata_objects.combine_objects.transpose().toList())

        // STEP 5: Collect versions
        versions = CISTOPIC_INFERCONSENSUS.out.versions
            .mix(
                CISTOPIC_QUALITYCONTROL.out.versions.first(),
                CISTOPIC_CREATEOBJECT.out.versions.first(),
                CISTOPIC_COMBINEOBJECTS.out.versions.ifEmpty{ Channel.empty() },
                ANNDATA_CONCAT.out.versions.ifEmpty{ Channel.empty() }
            )

        emit:
        consensus = consensus
        cistopic  = cistopic_objects.sample.mix(CISTOPIC_COMBINEOBJECTS.out.pkl)
        anndata   = anndata_objects.sample.mix(ANNDATA_CONCAT.out.h5ad)
        versions  = versions

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
    versions = Channel.empty()

    // Create pseudobulk, call peaks and update sample table
    if ( callPeaksFlag ) {
        PEAKCALLING(
            sample_table,
            celltypes,
            chromsizes,
            fragments_filename,
            narrowPeaks_dir
        )

        versions = versions.mix(PEAKCALLING.out.versions)
        // get pseudobulk peaks and updated sample table
        peak_metadata = PEAKCALLING.out.peaks
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
        versions = versions.mix(INFERPEAKS.out.versions)
    }

    emit:
    pseudobulk = callPeaksFlag ? PEAKCALLING.out.pseudobulk : Channel.empty()
    peaks      = callPeaksFlag ? PEAKCALLING.out.peaks : Channel.empty()
    consensus  = inferConsensusFlag ? INFERPEAKS.out.consensus : Channel.empty()
    cistopic   = inferConsensusFlag ? INFERPEAKS.out.cistopic : Channel.empty()
    anndata    = inferConsensusFlag ? INFERPEAKS.out.anndata : Channel.empty()
    versions   = versions
}