process CISTOPIC_QUALITYCONTROL {
    tag "Running QC for sample $meta.id"
    container 'docker://quay.io/cellgeni/pycistopic:2.0a0'
    
    input:
    tuple val(meta), path(fragments), path(fragments_index)
    tuple val(consensus_meta), path(consensus)
    tuple val(tss_bed_meta), path(tss_bed)

    output:
    tuple val(meta), path(fragments), path(fragments_index), path("qc/"), emit: qc
    path 'versions.yml',                                                  emit: versions
    
    script:
    def args = task.ext.args ?: ''
    """
    mkdir qc
    pycistopic qc \\
        $args \\
        --fragments $fragments \\
        --regions $consensus \\
        --tss $tss_bed \\
        --output qc/${meta.id} \\
        --threads $task.cpus
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pycisTopic: \$( python -c "import pycisTopic; print(pycisTopic.__version__)" )
        tss: ${tss_bed_meta.id}
    END_VERSIONS
    """
}