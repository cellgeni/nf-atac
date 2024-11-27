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
            --sample_table      specify a path to .csv file with sample names and path to the CellRanger-arc output dir (see example below)
            --celltypes         specify a path .csv file with celltype annotation or a path to pseudobulk_peaks.csv file with selected celltypes for consensus peak calling
            --callPeaks         if specified creates pseudobulks for celltypes specified in `--celltypes` for samples in `--sample_table`
            --inferConsensus    if specified runs a consensus peak calling and outputs cisTopic object for each sample in `--sample_table`

    Examples:
        1. Perform peak calling
            nextflow run main.nf --callPeaks --sample_table ./examples/sample_table.csv --celltypes ./examples/celltype_annotation.csv
        
        2. Infer consensus peaks and calculate features
            nextflow run main.nf --inferConsensus --sample_table ./results/updated_sample_table.csv --celltypes ./results/pseudobulk_peaks.csv
        
        3. Perform peak calling, infer consensus peaks and calculate features
            nextflow run main.nf --callPeaks --inferConsensus --sample_table ./examples/sample_table.csv --celltypes examples/celltype_annotation.csv

    == samples.csv format ==
    sample_id,outputdir
    WS_wEMB13386884,/lustre/path/to/cellranger-arc/output/
    WS_wEMB13386881,/lustre/path/to/cellranger-arc/output/

    == celltypes.csv format ==
    sample_id,barcode,celltype
    WS_wEMB13386884,AGAAGGTGTAATTAGC-1,vasculature
    WS_wEMB13386884,GATCGAGCACTTCATC-1,fibroblasts
    WS_wEMB13386881,ACAACATGTGATCAGC-1,vasculature
    WS_wEMB13386881,GAGCGGTCATGGAGGC-1,fibroblasts
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
    celltypes = file( params.celltypes )

    // Load other files required for cisTopic pipeline
    chromsizes = file( params.chromsizes )
    blacklist = file( params.blacklist )
    tss_bed = file( params.tss_bed )

    // Run PyCistopic pipeline
    PYCISTOPIC(
        sample_table,
        celltypes,
        chromsizes,
        blacklist,
        tss_bed,
        params.callPeaks,
        params.inferConsensus
    )

    // PYCISTOPIC.out.bed.view()
    // PYCISTOPIC.out.bigwig.view()
}

output {
    'narrowPeaks' {
        mode 'copy'
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