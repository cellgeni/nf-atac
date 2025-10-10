process ANNDATA_ATTACHCELLTYPES {
    tag "Attaching cell-type annotation to .obs section for ${meta.id}"
    container "quay.io/cellgeni/metacells-python:latest"
   
    input:
    tuple val(meta), path(h5ad, stageAs: 'input/*')
    path(metadata)
    
    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path '.attach_celltypes.log'   , emit: log
    path 'versions.yml'            , emit: versions
    
    script:
    barcode_column = task.ext.barcode_column ?: "obs_names"
    output = task.ext.output ?: "${meta.id}.h5ad"
    """
    attach_celltypes.py \
        --h5ad_file ${h5ad} \
        --sample_id ${meta.id} \
        --metadata ${metadata} \
        --barcode_column ${barcode_column} \
        --logfile .attach_celltypes.log \
        --output $output
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
    END_VERSIONS
    """
}