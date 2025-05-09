include { per_sample_preprocessing } from '../../modules/snapatac2'
include { combine_samples } from '../../modules/snapatac2'
include { call_peaks } from '../../modules/snapatac2'

workflow  SNAPATAC2 {
    take:
        sample_table

    main:
        samples = sample_table.splitCsv(skip: 1)
        
        per_sample_preprocessing(
          samples,
          params.min_counts,
          params.max_counts,
          params.min_tsse,
          params.n_features,
          params.genome
        )
        
        sample_file = per_sample_preprocessing.out.collectFile{item -> ["samples.csv", item[0]  + "," + item[1] + "," + item[2] + "\n"]}
        
        combine_samples(
          sample_file,
          params.n_features,
          params.genome
        )

        if(params.celltypes != null){
          call_peaks(combine_samples.out,
                     Channel.fromPath( params.celltypes ),
                     params.genome)
        }
}