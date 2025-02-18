// function to log error message if there is a memmory shortage
def lowMemoryError(sample, task_name) {
    log.warn "The memory is to low to perform ${task_name} for ${sample}"
    return 'retry'
}

// function to handle errors in MakePseudobulk process
def MakePseudobulkErrorHandler(exitStatus, celltypes) {
    if ( exitStatus == 130 ) {
        log.warn "The memory is to low to perform MakePseudobulk"
        return 'retry'
    } else if ( exitStatus == 0 ) {
        def celltype_name = celltypes.getName().split('\\.')[0]
        log.warn "No fragments found for ${celltype_name}"
        return 'ignore'
    } else {
        return 'terminate'
    }
}

// function to predict a memory consumption for MakePseudobulk process
def CalculatePseudobulkMemory(fragments_num, attempt) {
    def a = 2.3e-7
    def b = 32
    def x = fragments_num.toDouble()
    def mem = Math.ceil(a * x + b) + 32 * (attempt - 1)
    return mem
}

// Splits annotation file in separate files (per celltype) and calculates fragments numbers
process SplitCellTypeAnnotation {
    tag "Splitting celltype annotation"
    input:
        tuple val(sample_id_list), path(barcode_metrics, stageAs: 'barcode_metrics/*.csv')
        path(celltypes)
        path(sample_table)
    output:
        path('output/*.csv'), emit: celltypes
        path('fragments_celltype_x_sample.csv'), emit: fragments_celltype_x_sample
        path('fragments_per_celltype.csv'), emit: celltype_fragments
        path('updated_sample_table.csv'), emit: sample_table
        path('splitcelltypes.log'), emit: log
    script:
        def sample_id = sample_id_list.join(' ')
        """
        split_annotation.py \\
            --sample_id $sample_id \\
            --celltype_annotation $celltypes \\
            --barcode_metrics $barcode_metrics \\
            --sample_table $sample_table \\
            --output_dir output \\
            --fragments_celltype_x_sample fragments_celltype_x_sample.csv \\
            --fragments_per_celltype fragments_per_celltype.csv \\
            --updated_sample_table updated_sample_table.csv \\
            --logfile splitcelltypes.log \\
            --dropna
        """
    
}


// process to make pseudobulk profiles from fragments and celltype annototion
process MakePseudobulk {
    tag "Making pseudobulk for ${celltype_name}"
    input:
        tuple val(sample_id_list), path(fragments, stageAs: 'fragments/*/*'), path(fragments_index, stageAs: 'fragments/*/*')
        tuple val(celltype_name), path(celltypes), val(fragments_num)
        path(fragments_celltype_x_sample)
        path(chromsizes)
    output:
        tuple path('output/*.tsv.gz'), val(celltype_name), val(fragments_num), emit: pseudobulk_fragments
        path('bigwig/*.bw'), emit: bigwig
        path('*.log'), emit: log
    script:
        def sample_id = sample_id_list.join(' ')
        """
        make_pseudobulk.py \\
                --sample_id $sample_id \\
                --fragments $fragments \\
                --celltype_annotation $celltypes \\
                --fragments_celltype_x_sample $fragments_celltype_x_sample \\
                --chromsizes $chromsizes \\
                --output_dir output \\
                --bigwig_dir bigwig \\
                --cpus $task.cpus \\
                --logfile "pseudobulk.${celltype_name}.log"
        """

}

// Performs pseudobulk peak calling with MACS2
process PeakCalling {
    tag "Performing peak calling for ${celltype_name}"
    input:
        tuple path(pseudobulk_fragments), val(celltype_name), val(fragments_num)
        val(narrowPeaks_dir)
    output:
        tuple val(celltype_name), val(fragments_num), env(large_peaks_num), env(all_peaks_num), path("$narrowPeaks_dir/*.narrowPeak")
    script:
        """
        # perform peak calling
        peak_calling.py \\
                --bed_path $pseudobulk_fragments \\
                --output_dir $narrowPeaks_dir \\
                --cpus $task.cpus \\
                --skip_empty_peaks
        
        # count peaks
        large_peaks_num=\$(cut -f4 $narrowPeaks_dir/*.narrowPeak | sed 's/[a-z]\$//' | sort -u | wc -l)
        all_peaks_num=\$(wc -l < $narrowPeaks_dir/*.narrowPeak)
        """
}

// Make a csv table with 
process CollectPeakMetadata {
    input:
        tuple val(celltype_names), val(fragment_counts), val(large_peak_counts), val(all_peak_counts), path(narrowPeaks, stageAs: "$narrowPeaks_dir/*")
        val(narrowPeaks_dir)
    output:
        path('pseudobulk_peaks.tsv'), emit: metadata
        path(narrowPeaks), emit: narrow_peaks
    script:
        """
        collect_peak_info.py \\
                --publish_dir ${task.publishDir.path} \\
                --celltype_names ${celltype_names.join(' ')} \\
                --fragment_counts ${fragment_counts.join(' ')} \\
                --large_peak_counts ${large_peak_counts.join(' ')} \\
                --all_peak_counts ${all_peak_counts.join(' ')} \\
                --narrowPeaks ${narrowPeaks.join(' ')} \\
                --output pseudobulk_peaks.tsv
        """
}

// Derives consensus peaks from narrow peaks using the TGCA iterative peak filtering approach
process InferConsensus {
    tag "Infering consensus peaks"
    input:
        path("narrowPeaks/*")
        path(chromsizes)
        path(blacklist)
    output:
        path('*.bed'), emit: bed
        path('inferconsensus.log'), emit: log
    script:
        """
        infer_consensus.py \\
                --narrow_peaks ./narrowPeaks \\
                --chromsizes $chromsizes \\
                --blacklist $blacklist \\
                --consensus consensus_peaks.bed \\
                --skip_empty_peaks \\
                --logfile inferconsensus.log
        """
}

process QualityControl {
    tag "Running QC for sample $sample_id"
    input:
        tuple val(sample_id), path(fragments), path(fragments_index), val(fragments_num)
        path(consensus)
        path(tss_bed)
    output:
        tuple val(sample_id), path(fragments), path(fragments_index), path(consensus), path("qc/"), val(fragments_num)
    script:
        """
        mkdir qc
        pycistopic qc \\
            --fragments $fragments \\
            --regions $consensus \\
            --tss $tss_bed \\
            --output qc/$sample_id \\
            --threads $task.cpus
        """
}


process CreateCisTopicObject {
    tag "Creating cisTopic object for sample ${sample_id}"
    input:
        tuple val(sample_id), path(fragments), path(fragments_index), path(consensus), path(qc), val(fragments_num)
        path(blacklist)
    output:
        tuple val(sample_id),path("*.pkl")
        script:
        """
        create_cistopic.py \\
            --sample_id $sample_id \\
            --fragments $fragments \\
            --consensus $consensus \\
            --blacklist $blacklist \\
            --qc_dir $qc \\
            --cpus $task.cpus \\
            --use_automatic_thresholds
        """
}