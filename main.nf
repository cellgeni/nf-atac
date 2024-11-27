// NEXTFLOW FLAGS
nextflow.enable.dsl=2
nextflow.preview.output=true

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
    
    // Convert sample_table to path
    sample_table = file( params.sample_table )

    // Load celltype annotation file
    celltype_annotation = file( params.celltypes )

    // Load other files required for cisTopic pipeline
    chromsizes = file( params.chromsizes )
    blacklist = file( params.blacklist )
    tss_bed = file( params.tss_bed )

    // Run PyCistopic pipeline
    PYCISTOPIC(
        sample_table,
        celltype_annotation,
        chromsizes,
        blacklist,
        tss_bed,
        params.fromPseudobulk,
        params.cisTopicObject
    )

    // PYCISTOPIC.out.bed.view()
    // PYCISTOPIC.out.bigwig.view()
}

output {
    'narrowPeaks' {
        index {
            path '../pseudobulk_peaks.csv'
            header true
            mapper { celltype_name, fragments_num, large_peaks_num, all_peaks_num, narrow_peaks -> 
            [
                celltype: celltype_name,
                fragments: fragments_num,
                large_peaks: large_peaks_num,
                all_peaks: all_peaks_num,
                path: narrow_peaks
            ]
            }
        }
    }
}