// function to log error message if there is a memmory shortage
def lowMemoryError(sample, task_name) {
    log.warn "The memory is too low to perform ${task_name} for ${sample}"
    return 'retry'
}

process preprocess_sample {
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
  publishDir "${params.output_dir}/arrows_all/", mode: 'copy', pattern: '*.arrow'
  input:
    tuple val(sample_id), path(arrow_file)
    val(genome)
    val(n_samples)

  output:
    tuple path("gene_scores.mtx"), 
          path("gene_scores_obs.csv"), 
          path("gene_scores_var.csv"),
          val("gene_scores.h5ad"), emit: matrix_parts
    tuple val(sample_id), path(arrow_file), emit: arrows
  
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
    tuple path(mtx),path(obs),path(var),val(out_h5ad)
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

process call_peaks {
  publishDir "${params.output_dir}/arrows_sub/", mode: 'copy', pattern: '*.arrow'
  input:
    tuple val(sample_id), path(arrow_file, stageAs: 'input_arrows/*')
    path(celltype_file)
    val(genome)
    val(n_samples)

  output:
    tuple path("peak2cell.mtx"),
          path("peak2cell_obs.csv"),
          path("peak2cell_var.csv"),
          val("peak2cell.h5ad"), emit: matrix_parts
    path('*.arrow')
  
  script:
  """
  call_peaks.R \
   --arrow ${arrow_file.join(' ')} \
   --celltype_csv ${celltype_file} \
   --out_base peak2cell \
   --nthreads 8 \
   --genome $genome

   mv arrows/ArrowFiles/*arrow ./
  """
}