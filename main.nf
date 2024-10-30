// NEXTFLOW FLAGS
nextflow.enable.dsl = 2

// IMPORT SUBWORKFLOWS
include { PYCISTOPIC } from './workflows/pycistopic/main'

// HELP MESSAGE
def helpMessage() {
    log.info"""
    ===========================
    Peak Calling ATAC pipeline
    ===========================
    This pipeline performs peak calling for ATAC data (only cisTopic option is available at the moment)
    Usage: nextflow run main.nf [OPTIONS]
        options:
            --sample_table      specify a .csv file with sample names and path to the CellRanger-arc output dir (see example below)
            --celltypes         specify a .csv file with celltype annotation (must include 'sample_id', 'barcode' and 'celltype' colums)
            --cisTopicObject    if specified the script creates cisTopicObject for each sample (otherwise only `.bed` files with consensus peaks are available)

    Examples:
        1. Perform peak calling
            nextflow run main.nf --sample_table ./examples/samples.csv --celltypes examples/celltypes.csv
        
        2. Perform peak calling and additionally create a cisTopic object
            nextflow run main.nf --sample_table ./examples/samples.csv --celltypes examples/celltypes.csv --cisTopicObject

    == samples.csv format ==
    sample_id,outputdir
    WS_wEMB13386884,/lustre/path/to/cellranger-arc/output/
    WS_wEMB13386881,/lustre/path/to/cellranger-arc/output/
    ========================
    """.stripIndent()
}


// WORKFLOW
workflow {
    // Validate input arguments
    if (params.help) {
        helpMessage()
        System.exit(0)
    } else if (params.sample_table == null || params.celltypes == null) {
        helpMessage()
        error "Please specify all of the arguments listed above"
    }
    
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

    // Load celltype annotation file
    celltype_annotation = file( params.celltype_annotation )

    // Load other required files
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