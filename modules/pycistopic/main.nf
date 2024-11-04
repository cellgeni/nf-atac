// function to log error message if there is a memmory shortage
def lowMemoryError(sample, task_name) {
    log.warn "The memory is to low to perform ${task_name} for ${sample}"
    return 'retry'
}

process SplitCellTypeAnnotation {
    tag "Splitting celltype annotation"
    executor "local"
    input:
        path(sample_table)
        path(celltypes)
    output:
        path('*.csv')
    script:
        """
        # Get sample list
        tail -n +2 "${sample_table}" | cut -d "," -f1 | sort | uniq > samples.txt
        # Extract the header
        header=\$(head -n 1 "${celltypes}")
        # Get only relevant rows
        grep -f samples.txt $celltypes > filtered_annotation.csv

        # Read unique cell types, skipping the header
        tail -n +2 filtered_annotation.csv | cut -d "," -f3 | sort | uniq | while read celltype; do
            # Replace spaces with underscores in celltype for file naming
            safe_celltype=\$(echo "\${celltype}" | tr ' ' '_')
            
            # Create a file with the header for each cell type
            echo "\${header}" > "\${safe_celltype}.csv"
            
            # Append rows matching the cell type
            grep ",\$celltype\\\$" filtered_annotation.csv >> "\${safe_celltype}.csv"
        done
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
        tuple val(sample_id_list), path('output/*.tsv.gz'), emit: fragments
        path('output/*.bw'), emit: bigwig
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
                --skip_empty_fragments \\
                --cpus $task.cpus
        """

}

// Performs pseudobulk peak calling with MACS2
process PeakCalling {
    tag "Performing pseudobulk peak calling for sample $sample_id"
    input:
        tuple val(sample_id), path("fragments/*")
    output:
        tuple val(sample_id), path('narrowPeaks/*.narrowPeak')
    script:
        """
        peak_calling.py \\
                --bed_path ./fragments/ \\
                --output_dir narrowPeaks \\
                --cpus $task.cpus \\
                --skip_empty_peaks
        """
}

// Derives consensus peaks from narrow peaks using the TGCA iterative peak filtering approach
process InferConsensus {
    tag "Infering consensus peaks for sample $sample_id"
    input:
        tuple val(sample_id), path("narrowPeaks/*")
        path(chromsizes)
        path(blacklist)
    output:
        tuple val(sample_id), path('*.bed')
    script:
        """
        infer_consensus.py \\
                --sample_id $sample_id \\
                --narrow_peaks ./narrowPeaks \\
                --chromsizes $chromsizes \\
                --blacklist $blacklist \\
                --skip_empty_peaks
        """
}

process QualityControl {
    tag "Running QC for sample $sample_id"
    input:
        tuple val(sample_id), path(framgents), path(fragments_index), path(barcode_metrics), path(consensus)
        path(tss_bed)
    output:
        tuple val(sample_id), path(framgents), path(fragments_index), path(consensus), path("qc/")
    script:
        """
        mkdir qc
        pycistopic run qc \\
            --fragments $framgents \\
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