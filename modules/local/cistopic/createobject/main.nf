process CISTOPIC_CREATEOBJECT {
    tag "Creating cisTopic object for sample ${meta.id}"
    container '/nfs/cellgeni/singularity/images/scenicplus-fa55dae.sif'

    input:
    tuple val(meta), path(fragments), path(fragments_index), path(qc)
    tuple val(consensus_meta), path(consensus)
    tuple val(blacklist_meta), path(blacklist)
    
    output:
    tuple val(meta), path("*.pkl"), emit: pkl
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.json"), emit: json
    path 'versions.yml', emit: versions

    script:
    args = task.ext.args ?: ''
    """
    create_cistopic.py \\
        $args \\
        --sample_id ${meta.id} \\
        --fragments $fragments \\
        --consensus $consensus \\
        --blacklist $blacklist \\
        --qc_dir $qc \\
        --cpus $task.cpus \\
        --use_automatic_thresholds
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pycisTopic: \$( python -c "import pycisTopic; print(pycisTopic.__version__)" )
        blacklist: ${blacklist_meta.id}
    END_VERSIONS
    """
}