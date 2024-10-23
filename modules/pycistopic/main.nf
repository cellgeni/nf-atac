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
        path('output/*.tsv.gz'), emit: bed
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

process testProcess {
    output:
        path('*.csv')
    script:
        """
        make_pseudobulk.py \\
                --sample_id $sample_id
                --fragments $fragments
                --celltype-annotation $celltypes
                --chromsizes $chromsizes
                --cpus $task.cpus
        """
}

