include { preprocess_sample } from '../../modules/archr'
include { get_gene_scores } from '../../modules/archr'
include { make_h5ad } from '../../modules/archr'
include { call_peaks } from '../../modules/archr'

workflow  ARCHR {
    take:
        sample_table

    main:
        samples = sample_table.splitCsv(skip: 1).map{ sample_id, cellranger_output -> 
                                           [
                                                sample_id,
                                                file( "${cellranger_output}/${params.fragments_filename}" ),
                                                file( "${cellranger_output}/${params.fragments_filename}.tbi" )
                                           ]
                                           }
        
        preprocess_sample(
          samples,
          params.min_counts,
          params.min_tsse,
          params.genome,
        )

        arrows = preprocess_sample.out.map{sample_id, arrow_folder -> 
                  [
                    sample_id,
                    file("${arrow_folder}/${sample_id}.arrow")
                  ]}.toList().transpose().toList()
        
        gene_scores = get_gene_scores(
          arrows,
          params.genome,
          samples.count()
        )
        
        mtx_ch = gene_scores.matrix_parts

        if(params.celltypes != null){
          peaks = call_peaks(gene_scores.arrows,
                     Channel.fromPath( params.celltypes ),
                     params.genome,
                     samples.count())
          mtx_ch = mtx_ch.concat(peaks.matrix_parts)
        }
        
        mtx_ch.view()


        make_h5ad(mtx_ch,
                  samples.count())
                  
                  
}
