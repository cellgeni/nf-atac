// function to log error message if there is a memmory shortage
def lowMemoryError(sample, task_name) {
    log.warn "The memory is too low to perform ${task_name} for ${sample}"
    return 'retry'
}

process per_sample_preprocessing {
  publishDir "${params.output_dir}/", mode: 'copy'
  
  input:
    tuple val(sample_id), path(fragments)
    val(min_counts)
    val(max_counts)
    val(min_tsse)
    val(n_features)
    val(genome)
  
  
  output:
    tuple val(sample_id), path("${sample_id}.h5ad") //, path("${sample_id}_gene_matrix.h5ad")
  
  script:
  """
  preprocess_sample.py \
   --fragments $fragments \
   --sample_id $sample_id \
   --min_counts $min_counts \
   --max_counts $max_counts \
   --min_tsse $min_tsse \
   --n_features $n_features \
   --genome $genome
  """
}

process combine_samples {
  publishDir "${params.output_dir}/", mode: 'copy'
  
  input:
    path(samples)
    val(n_features)
    val(genome)

  output:
    tuple path("combined.h5ad"),path("combined_gene_matrix.h5ad")
  
  script:
  """
  combine_samples.py \
   --samples_file $samples \
   --n_features $n_features \
   --genome $genome \
   --postprocess
  """
}

process call_peaks {
  publishDir "${params.output_dir}/", mode: 'copy'
  
  input:
    tuple path(h5ad_file),path(gene_matrix)
    path(celltype_file)
    val(genome)

  output:
    path("*.h5ad")
    path('subadatas')
    path('merged_peaks.csv')
  
  script:
  """
  call_peaks.py \
   --h5ad_file $h5ad_file \
   --celltype_file $celltype_file \
   --genome $genome \
   --n_jobs 10
  """
}