process CISTOPIC_INFERCONSENSUS {
    tag "Infering consensus peaks"
    container '/nfs/cellgeni/singularity/images/scenicplus-fa55dae.sif'
    
    input:
    tuple val(meta), path("narrowPeaks/*")
    tuple val(meta_chromsizes), path(chromsizes)
    tuple val(meta_blacklist), path(blacklist)
   
    output:
    tuple val(meta), path('*.bed'), emit: bed
    path 'inferconsensus.log',      emit: log
    path 'versions.yml',            emit: versions

    
    script:
    def args = task.ext.args ?: ''
    """
    infer_consensus.py \\
        $args \\
        --narrow_peaks ./narrowPeaks \\
        --chromsizes $chromsizes \\
        --blacklist $blacklist \\
        --consensus consensus_peaks.bed \\
        --skip_empty_peaks \\
        --logfile inferconsensus.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pycisTopic: \$( python -c "import pycisTopic; print(pycisTopic.__version__)" )
        chromsizes: ${meta_chromsizes.id}
        blacklist: ${meta_blacklist.id}
    END_VERSIONS
    """
}