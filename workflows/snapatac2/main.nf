include { preprocess_sample } from '../../modules/snapatac2'
include { combine_samples } from '../../modules/snapatac2'
include { call_peaks } from '../../modules/snapatac2'

workflow  SNAPATAC2 {
    take:
        sample_table

    main:
        samples = sample_table.splitCsv(skip: 1)
        
        preprocess_sample(
          samples,
          params.min_counts,
          params.max_counts,
          params.min_tsse,
          params.n_features,
          params.genome,
          params.remove_doublets
        )
        
        combine_samples(
          preprocess_sample.out.toList().transpose().toList(),
          params.n_features,
          params.genome
        )

        if(params.celltypes != null){
          call_peaks(combine_samples.out,
                     Channel.fromPath( params.celltypes ),
                     params.genome)
        }
}