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

process CISTOPIC_PSEUDOBULK {
    tag "Making pseudobulk for ${celltype_meta.id}"
    container "docker://quay.io/cellgeni/pycistopic:2.0a0"
    
    input:
    tuple val(metalist), path(fragments, stageAs: 'sample_fragments/*/*'), path(fragments_index, stageAs: 'sample_fragments/*/*')
    tuple val(celltype_meta), path(celltypes)
    tuple val(chromsizes_meta), path(chromsizes)
    path(fragments_celltype_x_sample)
    
    output:
    tuple val(celltype_meta), path('fragments/*.tsv.gz'), emit: tsv
    path 'bigwig/*.bw',  emit: bigwig
    path '*.log',        emit: log
    path "versions.yml", emit: versions
    
    script:
    def samples = metalist.collect{ it.id }.join(' ')
    def args = task.ext.args ?: ''
    """
    make_pseudobulk.py \\
            $args \\
            --sample_id $samples \\
            --fragments $fragments \\
            --celltype_annotation $celltypes \\
            --fragments_celltype_x_sample $fragments_celltype_x_sample \\
            --chromsizes $chromsizes \\
            --output_dir fragments \\
            --bigwig_dir bigwig \\
            --cpus $task.cpus \\
            --logfile "pseudobulk.${celltype_meta.id}.log"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pycisTopic: \$( python -c "import pycisTopic; print(pycisTopic.__version__)" )
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
        chromsizes: ${chromsizes_meta.id}
    END_VERSIONS
    """
}