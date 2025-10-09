// IMPORT SUBWORKFLOWS
include { PYCISTOPIC } from './workflows/pycistopic/main'

// HELP MESSAGE
def helpMessage() {
    log.info"""
    ===========================
    Peak Calling ATAC pipeline
    ===========================
    This pipeline performs peak calling for ATAC data using pyCisTopic
    
    Usage: nextflow run main.nf [OPTIONS]
    
    Required arguments:
        --sample_table      Path to .csv file with sample names and paths to CellRanger-arc output directories
        --celltypes         Path to .csv file with celltype annotation OR path to pseudobulk_peaks.tsv file with selected celltypes
        --callPeaks         Run peak calling for provided celltypes (creates pseudobulks and calls peaks)
        --inferConsensus    Run consensus peak calling and feature calculation (creates cisTopic objects)
    
    Optional arguments:
        --output_dir        Output directory (default: 'results')
        --help              Show this help message
    
    Examples:
        1. Perform peak calling only:
            nextflow run main.nf --callPeaks --sample_table ./example/sample_table.csv --celltypes ./example/celltypes.csv
        
        2. Infer consensus peaks and calculate features (requires updated sample table from step 1):
            nextflow run main.nf --inferConsensus --sample_table ./example/updated_sample_table.csv --celltypes ./example/pseudobulk_peaks.tsv
        
        3. Run complete pipeline (peak calling + consensus inference):
            nextflow run main.nf --callPeaks --inferConsensus --sample_table ./example/sample_table.csv --celltypes ./example/celltypes.csv
    
    Input file formats:
        
        sample_table.csv:
        sample_id,outputdir
        WS_wEMB13386884,/path/to/cellranger-arc/output/
        WS_wEMB13386881,/path/to/cellranger-arc/output/
        
        celltypes.csv:
        sample_id,barcode,celltype
        WS_wEMB13386884,AGAAGGTGTAATTAGC-1,vasculature
        WS_wEMB13386884,GATCGAGCACTTCATC-1,fibroblasts
        WS_wEMB13386881,ACAACATGTGATCAGC-1,vasculature
        WS_wEMB13386881,GAGCGGTCATGGAGGC-1,fibroblasts
    
    For more information, see: https://github.com/cellgeni/nf-atac
    ========================
    """.stripIndent()
}

def lowResourceError(task_name) {
    log.warn "Not enough resources to perform ${task_name}"
    return 'retry'
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
    
    // Load files
    sample_table     = params.sample_table ? Channel.value( file( params.sample_table, checkIfExists: true ) ): Channel.empty()
    celltypes        = params.celltypes ? Channel.value( file( params.celltypes, checkIfExists: true ) ): Channel.empty()
    pseudobulk_peaks = params.pseudobulk_peaks ? Channel.value( file( params.pseudobulk_peaks, checkIfExists: true ) ): Channel.empty()
    atac_adata       = params.atac_adata ? Channel.value( file( params.atac_adata, checkIfExists: true ) ): Channel.empty()

    // Load other files required for cisTopic pipeline
    chromsizes = Channel.value( tuple( [id: "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"], file( params.chromsizes ) ) )
    blacklist  = Channel.value( tuple( [id: 'https://www.nature.com/articles/s41598-019-45839-z'], file( params.blacklist ) ) )
    tss_bed    = Channel.value( tuple( [id: 'https://github.com/cellgeni/nf-atac/blob/main/reference/hg38_pycistopic_tss.bed'], file( params.tss_bed ) ) )

    // Run PyCistopic pipeline
    PYCISTOPIC(
        sample_table,
        celltypes,
        pseudobulk_peaks,
        atac_adata,
        chromsizes,
        blacklist,
        tss_bed,
        params.callPeaks,
        params.inferConsensus,
        params.attachGEX
    )

    // Collect versions
    PYCISTOPIC.out.versions
        .splitText(by: 20)
        .unique()
        .collectFile(name: 'versions.yml', storeDir: params.output_dir, sort: true)
        .subscribe { __ -> 
            log.info("Versions saved to ${params.output_dir}/versions.yml")
        }
}