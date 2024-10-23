// process to make pseudobulk profiles from fragments and celltype annototion
process makePseudobulk {
    tag "Making pseudobulks for $sample_id"
    time { task.time + 30.minute * task.attempt }
    memory { task.memory + 50.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
    input:
        tuple val(sample_id), path(fragments), path(fragments_index)
        path(celltypes)
        path(chromsizes)
    output:
        tuple val(sample_id), path('output/*.tsv.gz'), emit: fragments
        path('output/*.bw'), emit: bigwig
    script:
        """
        make_pseudobulk.py \\
                --sample_id $sample_id \\
                --fragments $fragments \\
                --celltype_annotation $celltypes \\
                --chromsizes $chromsizes \\
                --cpus $task.cpus
        """

}

// Performs pseudobulk peak calling with MACS2
process peakCalling {
    tag "Performing pseudobulk peak calling for sample $sample_id"
    input:
        tuple val(sample_id), path("fragments/*")
    output:
        tuple val(sample_id), path('output/*.narrowPeak')
    script:
        """
        peak_calling.py \\
                --bed_path ./fragments/ \\
                --output_dir output \\
                --cpus $task.cpus
        """
}

// Derives consensus peaks from narrow peaks using the TGCA iterative peak filtering approach
process inferConsensus {
    debug true
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
                --blacklist $blacklist
        """
}