// function to log error message if there is a memmory shortage
def lowMemoryError(sample, task_name) {
    log.warn "The memory is too low to perform ${task_name} for ${sample}"
    return 'retry'
}

process preprocess_sample {
  input:
    tuple val(sample_id), path(fragments)
    val(min_counts)
    val(max_counts)
    val(min_tsse)
    val(n_features)
    val(genome)
    val(remove_doublets)
  
  
  output:
    tuple val(sample_id), path("${sample_id}.h5ad") 
  
  script:
  """
  preprocess_sample.py \
   --fragments $fragments \
   --sample_id $sample_id \
   --min_counts $min_counts \
   --max_counts $max_counts \
   --min_tsse $min_tsse \
   --n_features $n_features \
   --genome $genome \
   --remove_doublets $remove_doublets
  """
}

process combine_samples {
  publishDir "${params.output_dir}/full_adatas/", mode: 'copy'
  
  input:
    tuple val(sample_id), path(h5ad_in, stageAs: "anndatas/*")
    val(n_features)
    val(genome)

  output:
    tuple path("full.h5ads"),path("gene_matrix.h5ad"),path("anndatas")
  
  script:
  """
  combine_samples.py \
   --sample_id ${sample_id.join(' ')} \
   --h5ad_in $h5ad_in \
   --n_features $n_features \
   --genome $genome \
   --postprocess
  """
}

process call_peaks {
  publishDir "${params.output_dir}/", mode: 'copy'
  
  input:
    tuple path(h5ads),path(gene_matrix),path(anndatas, stageAs: "anndatas")
    path(celltype_file)
    val(genome)

  output:
    path('subset_adatas')
  
  script:
  """
  call_peaks.py \
   --h5ad_file $h5ads \
   --celltype_file $celltype_file \
   --genome $genome \
   --n_jobs 20
  """
}
