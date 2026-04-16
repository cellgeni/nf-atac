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
        --consensus         Path to consensus_peaks.bed (required when using --countPeaks without --inferConsensus)
        --atac_adata        Path to atac_anndata.csv (required when using --attachGEX without --countPeaks)
        
        :PIPELINE FILES:
        --chromsizes       Path to chromsizes file (default: reference/hg38.chrom.sizes)
        --blacklist        Path to blacklist file (default: reference/hg38-blacklist.v2.bed)
        --tss_bed          Path to TSS bed file (default: reference/hg38_pycistopic_tss.bed)

        :STEPS:
        --callPeaks         Run peak calling for provided celltypes (creates pseudobulks and calls peaks)
        --inferConsensus    Run consensus peak calling (creates consensus peaks)
        --countPeaks        Run QC and feature/object generation from consensus peaks (creates cisTopic and AnnData objects)
        --attachGEX         Attach GEX data to ATAC data for multiome integration

        :STEP RULES:
        --inferConsensus auto-enables --countPeaks when --countPeaks is not provided
        --countPeaks without --inferConsensus requires --consensus
        --attachGEX without --countPeaks requires --atac_adata
    
    Optional arguments:
        --output_dir        Output directory (default: 'results')
        --help              Show this help message
    
    Examples:
        1. Perform peak calling only:
            nextflow run main.nf --callPeaks --sample_table ./example/sample_table.csv --celltypes ./example/celltypes.csv
        
        2. Infer consensus peaks and calculate features (requires updated sample table from step 1):
            nextflow run main.nf --inferConsensus --countPeaks --sample_table ./results/updated_sample_table.csv --celltypes ./example/celltypes.csv --pseudobulk_peaks ./results/pseudobulk_peaks.csv

        3. Run countPeaks from an existing consensus file:
            nextflow run main.nf --countPeaks --sample_table ./results/updated_sample_table.csv --celltypes ./example/celltypes.csv --consensus ./results/consensus_peaks.bed
        
        4. Attach GEX data to existing ATAC data:
            nextflow run main.nf --attachGEX --sample_table ./example/updated_sample_table.csv --celltypes ./example/celltypes.csv --atac_adata ./results/atac_anndata.csv
        
        5. Infer consensus peaks, count peaks, and attach GEX in one go:
            nextflow run main.nf --inferConsensus --countPeaks --attachGEX --sample_table ./results/updated_sample_table.csv --celltypes ./example/celltypes.csv --pseudobulk_peaks ./results/pseudobulk_peaks.csv

        6. Run complete pipeline (peak calling + consensus inference + countPeaks + GEX attachment):
            nextflow run main.nf --callPeaks --inferConsensus --countPeaks --attachGEX --sample_table ./example/sample_table.csv --celltypes ./example/celltypes.csv
    
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

// WORKFLOW
workflow {
    // Validate input arguments
    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    // Normalize step flags.
    callPeaksFlag       = params.callPeaks
    inferConsensusFlag  = params.inferConsensus
    countPeaksFlag      = params.countPeaks
    attachGEXFlag       = params.attachGEX

    // Check required arguments for peak calling
    if ( callPeaksFlag && (! params.sample_table || ! params.celltypes || ! params.chromsizes) ) {
        error("Please provide --sample_table, --celltypes and --chromsizes when using --callPeaks")
    }

    // Check required arguments for consensus peak inference
    if ( inferConsensusFlag && (! params.chromsizes || ! params.blacklist) ) {
        error("Please provide --chromsizes and --blacklist when using --inferConsensus")
    }

    // Check required arguments for countPeaks
    if ( countPeaksFlag && (! params.sample_table || ! params.celltypes || ! params.blacklist || ! params.tss_bed) ) {
        error("Please provide --sample_table, --celltypes, --blacklist and --tss_bed when using --countPeaks")
    }

    // countPeaks requires consensus peaks from either inferConsensus step or external file
    if ( countPeaksFlag && ! inferConsensusFlag && ! params.consensus ) {
        error("Please provide --consensus when using --countPeaks without --inferConsensus")
    }

    // Check required arguments for attaching GEX data
    if ( attachGEXFlag && (! params.sample_table || ! params.celltypes) ) {
        error("Please provide --sample_table and --celltypes when using --attachGEX")
    }

    // Check required arguments for INFERPEAKS without CALLPEAKS
    if ( inferConsensusFlag && ! callPeaksFlag && ! params.pseudobulk_peaks ) {
        error("Please provide --pseudobulk_peaks when using --inferConsensus without --callPeaks")
    }

    // Validate pseudobulk_peaks CSV has required columns: celltype, fragments, path
    if ( inferConsensusFlag && ! callPeaksFlag && params.pseudobulk_peaks ) {
        def pbp_file = file( params.pseudobulk_peaks )
        def header = pbp_file.readLines().first().split(',').collect { col -> col.trim() }
        def required_cols = ['celltype', 'path']
        def missing_cols = required_cols.findAll { col -> ! header.contains(col) }
        if ( missing_cols ) {
            error("--pseudobulk_peaks CSV is missing required columns: ${missing_cols.join(', ')}. Expected columns: ${required_cols.join(', ')}")
        }
    }

    // AttachGEX needs ATAC anndata from countPeaks output or external atac_adata file
    if ( attachGEXFlag && ! countPeaksFlag && ! params.atac_adata ) {
        error("Please provide --atac_adata when using --attachGEX without --countPeaks")
    }

    // Load files
    sample_table     = params.sample_table ? channel.value( file( params.sample_table, checkIfExists: true ) ): channel.empty()
    celltypes        = params.celltypes ? channel.value( file( params.celltypes, checkIfExists: true ) ): channel.empty()
    pseudobulk_peaks = params.pseudobulk_peaks ? channel.value( file( params.pseudobulk_peaks, checkIfExists: true ) ): channel.empty()
    consensus        = params.consensus ? channel.value( [ [id: "input_consensus"], file( params.consensus, checkIfExists: true ) ] ): channel.empty()
    atac_adata       = params.atac_adata ? channel.value( file( params.atac_adata, checkIfExists: true ) ): channel.empty()

    // Load other files required for cisTopic pipeline
    chromsizes = channel.value( tuple( [id: "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"], file( params.chromsizes ) ) )
    blacklist  = channel.value( tuple( [id: 'https://www.nature.com/articles/s41598-019-45839-z'], file( params.blacklist ) ) )
    tss_bed    = channel.value( tuple( [id: 'https://github.com/cellgeni/nf-atac/blob/main/reference/hg38_pycistopic_tss.bed'], file( params.tss_bed ) ) )

    // Run PyCistopic pipeline
    PYCISTOPIC(
        sample_table,
        celltypes,
        pseudobulk_peaks,
        consensus,
        atac_adata,
        chromsizes,
        blacklist,
        tss_bed,
        callPeaksFlag,
        inferConsensusFlag,
        countPeaksFlag,
        attachGEXFlag,
        params.cistopic.gex_filtered
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
    def versions_header = "\"NF-ATAC\":\n    name: ${workflow.manifest.name}\n    version: ${workflow.manifest.version}\n    description: \"${workflow.manifest.description}\"\n"
    PYCISTOPIC.out.versions
        .splitText(by: 20)
        .unique()
        .collectFile(name: 'versions.yml', storeDir: params.output_dir, sort: true, seed: versions_header)
        .subscribe { __ -> 
            log.info("Versions saved to ${params.output_dir}/versions.yml")
        }
}