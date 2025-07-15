// function to log error message if there is a memmory shortage
def lowMemoryError(sample, task_name) {
    log.warn "The memory is too low to perform ${task_name} for ${sample}"
    return 'retry'
}

process preprocess_sample {
  publishDir "${params.output_dir}/arrows/", mode: 'copy'
  input:
    tuple val(sample_id), path(fragments), path(fragments_index)
    val(min_counts)
    val(min_tsse)
    val(genome)

  
  
  output:
    tuple val(sample_id), path("${sample_id}"), emit: arrows
  
  script:
  """
  create_arrow_file.R \
   --fragments $fragments \
   --sample_id $sample_id \
   --min_counts $min_counts \
   --min_tsse $min_tsse \
   --genome $genome \
   --nthreads 8
  """
}

process get_gene_scores {
  input:
    tuple val(sample_id), path(arrow_file)
    val(genome)
    val(n_samples)

  output:
    path("gene_scores.mtx"), emit: mtx
    path("gene_scores_obs.csv"), emit: obs
    path("gene_scores_var.csv"), emit: var
  
  script:
  """
  get_gene_scores.R \
   --arrow ${arrow_file.join(' ')} \
   --out_base gene_scores \
   --nthreads 8 \
   --genome $genome
  """
}

process make_h5ad{
  publishDir "${params.output_dir}/", mode: 'copy'
  input:
    path(mtx)
    path(obs)
    path(var)
    val(out_h5ad)
    val(n_samples)

  output: 
    path(out_h5ad), emit: h5ad

  script:
  """
  make_h5ad.py \
    --mtx ${mtx} \
    --obs ${obs} \
    --var ${var} \
    --out_h5ad ${out_h5ad}
  """
}

/*process combine_samples {
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
    path(blacklist)

  output:
    path('subset_adatas')
  
  script:
  """
  call_peaks.py \
   --h5ad_file $h5ads \
   --celltype_file $celltype_file \
   --genome $genome \
   --blacklist $blacklist \
   --n_jobs 20
  """
}
*/