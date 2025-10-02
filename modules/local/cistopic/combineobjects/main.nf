process CISTOPIC_COMBINEOBJECTS {
    tag "Combining Python Objects"
    container '/nfs/cellgeni/singularity/images/scenicplus-fa55dae.sif'

    input:
    tuple val(metalist), path(cistopic, stageAs: "cistopic/*.pkl")
    
    output:
    tuple val(metalist), path("combined_cistopic_object.pkl"), emit: pkl
    path 'versions.yml', emit: versions
    
    script:
    """
    combine_objects.py \\
        --cistopic $cistopic \\
        --combined_cistopic combined_cistopic_object.pkl
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pycisTopic: \$( python -c "import pycisTopic; print(pycisTopic.__version__)" )
    END_VERSIONS
    """

}