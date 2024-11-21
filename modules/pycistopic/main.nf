// function to log error message if there is a memmory shortage
def lowMemoryError(sample, task_name) {
    log.warn "The memory is to low to perform ${task_name} for ${sample}"
    return 'retry'
}

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

def MakePseudobulkMemory(fragments_num, attempt) {
    def a = 2.3e-7
    def b = 32
    def x = fragments_num.toDouble()
    def mem = Math.ceil(a * x + b) + 32 * (attempt - 1)
    return "${mem as int}.GB"
}

process SplitCellTypeAnnotation {
    tag "Splitting celltype annotation"
    input:
        tuple val(sample_id_list), path(fragments, stageAs: 'fragments/*/*'), path(fragments_index, stageAs: 'fragments/*/*'), path(barcode_metrics, stageAs: 'fragments/*/*')
        path(celltypes)
    output:
        path('output/*.csv'), emit: celltypes
        path('fragments_celltype_x_sample.csv'), emit: fragments_celltype_x_sample
        path('fragments_per_celltype.csv'), emit: celltype_fragments
        path('fragments_per_sample.csv'), emit: sample_fragments
        path('splitcelltypes.log'), emit: log
    script:
        def sample_id = sample_id_list.join(' ')
        """
        split_annotation.py \\
            --sample_id $sample_id \\
            --celltype_annotation $celltypes \\
            --barcode_metrics $barcode_metrics \\
            --output_dir output \\
            --filtered_sample_table filtered_sample_table.csv \\
            --fragments_celltype_x_sample fragments_celltype_x_sample.csv \\
            --fragments_per_celltype fragments_per_celltype.csv \\
            --fragments_per_sample fragments_per_sample.csv \\
            --logfile splitcelltypes.log \\
            --dropna
        """
    
}


// process to make pseudobulk profiles from fragments and celltype annototion
process MakePseudobulk {
    tag "Making pseudobulk for ${celltypes.getName().split('\\.')[0]}"
    input:
        tuple val(sample_id_list), path(fragments, stageAs: 'fragments/*/*'), path(fragments_index, stageAs: 'fragments/*/*'), path(barcode_metrics, stageAs: 'fragments/*/*')
        tuple val(celltype_name), path(celltypes), val(fragments_num)
        path(chromsizes)
    output:
        path('output/*.tsv.gz'), emit: fragments
        path('bigwig/*.bw'), emit: bigwig
        path('*.log'), emit: log
    script:
        def sample_id = sample_id_list.join(' ')
        """
        make_pseudobulk.py \\
                --sample_id $sample_id \\
                --fragments $fragments \\
                --celltype_annotation $celltypes \\
                --chromsizes $chromsizes \\
                --barcode_metrics $barcode_metrics \\
                --output_dir output \\
                --bigwig_dir bigwig \\
                --skip_empty_fragments \\
                --cpus $task.cpus \\
                --logfile "pseudobulk.${celltype_name}.log"
        """

}

// Performs pseudobulk peak calling with MACS2
process PeakCalling {
    tag "Performing pseudobulk peak calling for sample ${fragments.getName()}"
    input:
        path(fragments)
    output:
        path('narrowPeaks/*.narrowPeak')
    script:
        """
        peak_calling.py \\
                --bed_path $fragments \\
                --output_dir narrowPeaks \\
                --cpus $task.cpus \\
                --skip_empty_peaks
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
        tuple val(sample_id), path(fragments), path(fragments_index), path(barcode_metrics)
        path(consensus)
        path(tss_bed)
    output:
        tuple val(sample_id), path(fragments), path(fragments_index), path(consensus), path("qc/")
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
        tuple val(sample_id), path(fragments), path(fragments_index), path(consensus), path(qc)
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