process ANNDATA_COUPLEMULTIOME {
    tag "Entangling ${meta.id}'s multiome data"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/metacells-python:latest':
        'quay.io/cellgeni/metacells-python:latest' }"

    input:
    tuple val(meta), path(atac), path(gex)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.h5mu"), emit: h5mu
    path '.couple_multiome.log'    , emit: log
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_coupled"
    """
    couple_multiome.py \
            --gex ${gex} \
            --atac ${atac} \
            --prefix ${prefix} \
            --logfile .couple_multiome.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        mudata: \$( python -c "import mudata; print(mudata.__version__)" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_coupled"
    """
    touch .couple_multiome.log
    touch ${prefix}_gex.h5ad
    touch ${prefix}_atac.h5ad
    touch ${prefix}.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        mudata: \$( python -c "import mudata; print(mudata.__version__)" )
    END_VERSIONS
    """
}