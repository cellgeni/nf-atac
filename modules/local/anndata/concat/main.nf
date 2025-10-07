process ANNDATA_CONCAT {
    tag "Combining Python Objects"
    container 'docker://quay.io/cellgeni/toh5ad:latest'

    input:
    tuple val(metalist), path(anndata, stageAs: "anndata/*.h5ad")
    
    output:
    tuple val(metalist), path("*.h5ad"), emit: h5ad
    path 'versions.yml', emit: versions
    
    script:
    def output = task.ext.output ?: 'combined_anndata_object.h5ad'
    """
    concat_adata.py \\
        --anndata $anndata \\
        --output $output
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
    END_VERSIONS
    """

}