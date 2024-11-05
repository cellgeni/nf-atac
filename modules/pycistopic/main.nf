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

process SplitCellTypeAnnotation {
    tag "Splitting celltype annotation"
    input:
        path(sample_table)
        path(celltypes)
    output:
        path('output/*.csv'), emit: celltypes
        path('filtered_sample_table.csv'), emit: sample_table
        path('pseudobulk.log'), emit: log
    script:
        """
        split_annotation.py \\
            --sample_table $sample_table \\
            --celltype_annotation $celltypes \\
            --output_dir output \\
            --filtered_sample_table filtered_sample_table.csv \\
            --logfile splitcelltypes.log \\
            --dropna
        """
    
}


// process to make pseudobulk profiles from fragments and celltype annototion
process MakePseudobulk {
    tag "Making pseudobulks for ${sample_id_list.size} samples"
    input:
        tuple val(sample_id_list), path(fragments, stageAs: 'fragments/*/*'), path(fragments_index, stageAs: 'fragments/*/*'), path(barcode_metrics, stageAs: 'fragments/*/*')
        each path(celltypes)
        path(chromsizes)
    output:
        path('output/*.tsv.gz'), emit: fragments
        path('bigwig/*.bw'), emit: bigwig
        path('pseudobulk.log'), emit: log
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
                --logfile pseudobulk.log
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
        pycistopic qc run \\
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