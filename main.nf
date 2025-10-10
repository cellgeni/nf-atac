// IMPORT SUBWORKFLOWS
include { PYCISTOPIC } from './workflows/pycistopic/main'

// HELP MESSAGE
def helpMessage() {
    log.info"""
    ===========================
    Peak Calling ATAC pipeline
    ===========================
    This pipeline performs peak calling for ATAC data using pyCisTopic and supports multiome data integration
    
    Usage: nextflow run main.nf [OPTIONS]
    
    Required arguments:
        :INPUT FILES:
        --sample_table      Path to .csv file with sample names and paths to CellRanger-arc output directories
        --celltypes         Path to .csv file with celltype annotation
        --pseudobulk_peaks  Path to pseudobulk_peaks.csv (required when using --inferConsensus without --callPeaks)
        --atac_adata        Path to atac_anndata.csv (required when using --attachGEX without --inferConsensus)
        
        :PIPELINE FILES:
        --chromsizes       Path to chromsizes file (default: reference/hg38.chrom.sizes)
        --blacklist        Path to blacklist file (default: reference/hg38-blacklist.v2.bed)
        --tss_bed          Path to TSS bed file (default: reference/hg38_pycistopic_tss.bed)

        :STEPS:
        --callPeaks         Run peak calling for provided celltypes (creates pseudobulks and calls peaks)
        --inferConsensus    Run consensus peak calling and feature calculation (creates cisTopic objects)
        --attachGEX         Attach GEX data to ATAC data for multiome integration
    
    Optional arguments:
        --output_dir        Output directory (default: 'results')
        --help              Show this help message
    
    Examples:
        1. Perform peak calling only:
            nextflow run main.nf --callPeaks --sample_table ./example/sample_table.csv --celltypes ./example/celltypes.csv
        
        2. Infer consensus peaks and calculate features (requires updated sample table from step 1):
            nextflow run main.nf --inferConsensus --sample_table ./results/updated_sample_table.csv --celltypes ./example/celltypes.csv --pseudobulk_peaks ./results/pseudobulk_peaks.csv
        
        3. Attach GEX data to existing ATAC data:
            nextflow run main.nf --attachGEX --sample_table ./example/updated_sample_table.csv --celltypes ./example/celltypes.csv --atac_adata ./results/atac_anndata.csv
        
        4. Infer consensus peaks and attach GEX in one go:
            nextflow run main.nf --inferConsensus --attachGEX --sample_table ./results/updated_sample_table.csv --celltypes ./example/celltypes.csv --pseudobulk_peaks ./results/pseudobulk_peaks.csv

        5. Run complete pipeline (peak calling + consensus inference + GEX attachment):
            nextflow run main.nf --callPeaks --inferConsensus --attachGEX --sample_table ./example/sample_table.csv --celltypes ./example/celltypes.csv
    
    Input file formats:
        
        sample_table.csv (basic format):
        sample_id,path
        WS_wEMB13386884,/path/to/cellranger-arc/output/
        WS_wEMB13386881,/path/to/cellranger-arc/output/
        
        updated_sample_table.csv (with fragment counts):
        sample_id,path,fragments
        WS_wEMB13386884,/path/to/cellranger-arc/output/,730872409
        WS_wEMB13386881,/path/to/cellranger-arc/output/,1118846819
        
        celltypes.csv:
        sample_id,barcode,celltype
        WS_wEMB13386884,AGAAGGTGTAATTAGC-1,vasculature
        WS_wEMB13386884,GATCGAGCACTTCATC-1,fibroblasts
        WS_wEMB13386881,ACAACATGTGATCAGC-1,vasculature
        WS_wEMB13386881,GAGCGGTCATGGAGGC-1,fibroblasts
        
        atac_anndata.csv:
        sample_id,path
        WS_wEMB13386884,/path/to/WS_wEMB13386884.h5ad
        WS_wEMB13386881,/path/to/WS_wEMB13386881.h5ad
    
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
    }

    // Check required arguments for peak calling
    if ( params.callPeaks && (! params.sample_table || ! params.celltypes || ! params.chromsizes) ) {
        error("Please provide --sample_table, --celltypes and --chromsizes when using --callPeaks")
    }

    // Check required arguments for consensus peak inference
    if ( params.inferConsensus && (! params.sample_table || ! params.celltypes || ! params.chromsizes || ! params.blacklist || ! params.tss_bed) ) {
        error("Please provide --sample_table, --celltypes, --chromsizes, --blacklist and --tss_bed when using --inferConsensus")
    }

    // Check required arguments for attaching GEX data
    if ( params.attachGEX && (! params.sample_table || ! params.celltypes) ) {
        error("Please provide --sample_table and --celltypes when using --attachGEX")
    }

    // Check required arguments for INFERPEAKS without CALLPEAKS
    if ( params.inferConsensus && ! params.callPeaks && ! params.pseudobulk_peaks ) {
        error("Please provide --pseudobulk_peaks when using --inferConsensus without --callPeaks")
    }

    // Check required arguments for ATTACHGEX without INFERPEAKS
    if ( params.attachGEX && ! params.inferConsensus && ! params.atac_adata ) {
        error("Please provide --atac_adata when using --attachGEX without --inferConsensus")
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

    // Collect ATAC anndata object paths (if generated)
    PYCISTOPIC.out.atac_anndata
        .collectFile(
            name: 'atac_anndata.csv',
            storeDir: params.output_dir,
            newLine: true,
            seed: "sample_id,path",
            sort: true
        ) { meta, path ->
            "${meta.id},${path.toString()}"
        }
        .subscribe { __ -> 
            log.info("ATAC anndata paths saved to ${params.output_dir}/atac_anndata.csv")
        }

    // Collect versions
    PYCISTOPIC.out.versions
        .splitText(by: 20)
        .unique()
        .collectFile(name: 'versions.yml', storeDir: params.output_dir, sort: true)
        .subscribe { __ -> 
            log.info("Versions saved to ${params.output_dir}/versions.yml")
        }
}