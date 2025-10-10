process ANNDATA_TOH5AD {
    tag "Converting ${meta.id}'s file to .h5ad"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/metacells-python:latest':
        'quay.io/cellgeni/metacells-python:latest' }"

    input:
    tuple val(meta), path(input, name: 'input/*')

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convert_to_h5ad.py \
            --input ${input} \
            --sample_id ${meta.id} \
            ${args} \
            --output ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
    END_VERSIONS
    """
}