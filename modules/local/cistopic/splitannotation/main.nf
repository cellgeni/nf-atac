// Splits annotation file in separate files (per celltype) and calculates fragments numbers
process CISTOPIC_SPLITANNOTATION {
    tag "Splitting celltype annotation"
    input:
    tuple val(metalist), path(barcode_metrics, stageAs: 'barcode_metrics/*.csv')
    path(celltypes)
    path(sample_table)
    
    output:
    path 'output/*.csv',                    emit: celltypes
    path 'fragments_celltype_x_sample.csv', emit: fragments_celltype_x_sample
    path 'fragments_per_celltype.csv',      emit: celltype_fragments
    path 'updated_sample_table.csv',        emit: sample_table
    path 'splitcelltypes.log',              emit: log
    path "versions.yml",                    emit: versions
    
    script:
    def samples = metalist.collect{ it.id }.join(' ')
    """
    split_annotation.py \\
        --sample_id $samples \\
        --celltype_annotation $celltypes \\
        --barcode_metrics $barcode_metrics \\
        --sample_table $sample_table \\
        --output_dir output \\
        --fragments_celltype_x_sample fragments_celltype_x_sample.csv \\
        --fragments_per_celltype fragments_per_celltype.csv \\
        --updated_sample_table updated_sample_table.csv \\
        --logfile splitcelltypes.log \\
        --dropna
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        scanpy: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """
    
}