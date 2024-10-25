// NEXTFLOW FLAGS
nextflow.enable.dsl = 2

// IMPORT SUBWORKFLOWS
include { PYCISTOPIC } from './workflows/pycistopic/main'

// WORKFLOW
workflow {
    // Read sample table and convert path to dir to path to fragments and barcode matrics
    sample_table = Channel
                        .fromPath(params.sample_table, checkIfExists: true)
                        .splitCsv(skip: 1)
                        .map{ sample_id, cellranger_arc_output -> 
                                [
                                    sample_id,
                                    file( "${cellranger_arc_output}/${params.fragments_filename}" ),
                                    file( "${cellranger_arc_output}/${params.fragments_filename}.tbi" ),
                                    file( "${cellranger_arc_output}/${params.barcode_metrics_filename}" ),
                                ]
                            }

    // Load input files
    celltype_annotation = file( params.celltype_annotation )
    chromsizes = file( params.chromsizes )
    blacklist = file( params.blacklist )
    tss_bed = file( params.tss_bed )

    // Run PyCistopic pipeline
    PYCISTOPIC(
        sample_table,
        celltype_annotation,
        chromsizes,
        blacklist,
        tss_bed
    )

    // PYCISTOPIC.out.bed.view()
    // PYCISTOPIC.out.bigwig.view()
}