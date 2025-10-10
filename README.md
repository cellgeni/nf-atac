# ATAC pipeline
This pipeline performs peak calling for ATAC-seq data using pyCisTopic. It supports pseudobulk creation, peak calling, consensus peak inference, cisTopic object generation, and multiome integration with GEX data for downstream analysis.

## Contents of Repo
* `main.nf` - the Nextflow pipeline that runs the whole pipeline
* `modules/local/` - a collection of pipeline processes organized by component
  * `modules/local/cistopic/` - cisTopic-related processes (countfragments, splitannotation, pseudobulk, callpeaks, inferconsensus, qualitycontrol, createobject, combineobjects)
  * `modules/local/anndata/` - AnnData-related processes (concat, attachcelltypes, couplemultiome, toh5ad)
* `workflows/` - high-level workflow orchestration
  * `workflows/pycistopic/main.nf` - main pyCisTopic workflow
* `subworkflows/local/` - reusable subworkflow components
  * `subworkflows/local/cistopic_peakcalling/` - peak calling subworkflow
  * `subworkflows/local/cistopic_inferpeaks/` - consensus peak inference subworkflow
  * `subworkflows/local/cistopic_attachgex/` - GEX data attachment subworkflow
* `configs/` - modular configuration files for each process
* `reference/` - reference files (chromosome sizes, blacklists, TSS annotations)
* `nextflow.config` - the main configuration script

## Pipeline Arguments

### Required Arguments

#### Input Files
* `--sample_table`: Path to .csv file with sample names and paths to the CellRanger-arc output directories
* `--celltypes`: Path to .csv file with celltype annotation
* `--pseudobulk_peaks`: Path to pseudobulk_peaks.csv (required when using --inferConsensus without --callPeaks)
* `--atac_adata`: Path to atac_anndata.csv (required when using --attachGEX without --inferConsensus)

#### Pipeline Files
* `--chromsizes`: Path to chromsizes file (default: reference/hg38.chrom.sizes)
* `--blacklist`: Path to blacklist file (default: reference/hg38-blacklist.v2.bed)
* `--tss_bed`: Path to TSS bed file (default: reference/hg38_pycistopic_tss.bed)

#### Steps
* `--callPeaks`: Run peak calling for provided celltypes (creates pseudobulks and calls peaks)
* `--inferConsensus`: Run consensus peak calling and feature calculation (creates cisTopic objects)
* `--attachGEX`: Attach GEX data to ATAC data for multiome integration

### Optional Arguments
* `--output_dir`: Output directory (default: 'results')
* `--help`: Show help message

## Pipeline Parameters

### Core Pipeline Parameters
| Parameter        | Default     | Description                                             |
| ---------------- | ----------- | ------------------------------------------------------- |
| `--output_dir`   | `'results'` | Output directory for all pipeline results               |
| `--publish_mode` | `'link'`    | How to publish output files ('link', 'copy', 'symlink') |

### Reference Files
| Parameter      | Default                             | Description                                           |
| -------------- | ----------------------------------- | ----------------------------------------------------- |
| `--chromsizes` | `reference/hg38.chrom.sizes`        | Chromosome sizes file for hg38 genome                 |
| `--blacklist`  | `reference/hg38-blacklist.v2.bed`   | ENCODE blacklist regions to exclude                   |
| `--tss_bed`    | `reference/hg38_pycistopic_tss.bed` | Transcription start sites for TSS enrichment analysis |

### pyCisTopic Parameters

#### Pseudobulking
| Parameter                     | Default | Description                                        |
| ----------------------------- | ------- | -------------------------------------------------- |
| `--cistopic.normalize_bigwig` | `true`  | Generate normalized BigWig files for visualization |

#### Peak Calling
| Parameter                     | Default   | Description                                               |
| ----------------------------- | --------- | --------------------------------------------------------- |
| `--cistopic.specie`           | `'hs'`    | Species for peak calling ('hs' for human, 'mm' for mouse) |
| `--cistopic.input_format`     | `'BEDPE'` | Input format for MACS2 peak calling                       |
| `--cistopic.shift`            | `73`      | Shift size for MACS2 peak calling                         |
| `--cistopic.extend_read_size` | `146`     | Fragment extension size for MACS2                         |
| `--cistopic.keep_duplicates`  | `'all'`   | How to handle duplicate reads ('all', 'auto', 'none')     |
| `--cistopic.q_value_cutoff`   | `0.05`    | Q-value threshold for peak calling                        |

#### Consensus Peak Inference
| Parameter                    | Default | Description                                  |
| ---------------------------- | ------- | -------------------------------------------- |
| `--cistopic.peak_half_width` | `250`   | Half-width for resizing consensus peaks (bp) |

#### Quality Control
| Parameter                                 | Default | Description                                            |
| ----------------------------------------- | ------- | ------------------------------------------------------ |
| `--cistopic.tss_flank_window`             | `2000`  | Window size around TSS for enrichment calculation (bp) |
| `--cistopic.tss_smoothing_rolling_window` | `10`    | Rolling window for TSS enrichment smoothing            |
| `--cistopic.tss_window`                   | `50`    | Central window size for TSS enrichment                 |
| `--cistopic.tss_min_norm`                 | `0.2`   | Minimum normalization value for TSS enrichment         |
| `--cistopic.min_fragments_per_cb`         | `10`    | Minimum fragments per cell barcode for QC              |
| `--cistopic.use_pyranges`                 | `false` | Use PyRanges for genomic operations                    |
| `--cistopic.dont_collapse_duplicates`     | `false` | Preserve duplicate fragments in analysis               |

#### cisTopic Object Creation
| Parameter                             | Default | Description                                            |
| ------------------------------------- | ------- | ------------------------------------------------------ |
| `--cistopic.min_frag`                 | `1`     | Minimum fragments per peak for inclusion               |
| `--cistopic.min_cell`                 | `1`     | Minimum cells per peak for inclusion                   |
| `--cistopic.is_acc`                   | `1`     | Whether data is accessibility data (1=yes, 0=no)       |
| `--cistopic.split_pattern`            | `'___'` | Pattern for splitting cell identifiers                 |
| `--cistopic.check_for_duplicates`     | `true`  | Check for duplicate peaks in consensus                 |
| `--cistopic.use_automatic_thresholds` | `true`  | Use automatic QC thresholds based on data distribution |

### Parameter Usage Examples
```bash
# Override default output directory
nextflow run main.nf --callPeaks --sample_table sample.csv --celltypes celltypes.csv --output_dir my_results

# Adjust peak calling parameters
nextflow run main.nf --callPeaks --sample_table sample.csv --celltypes celltypes.csv \
  --cistopic.q_value_cutoff 0.01 --cistopic.peak_half_width 300

# Use custom reference files
nextflow run main.nf --callPeaks --sample_table sample.csv --celltypes celltypes.csv \
  --chromsizes /path/to/custom.chrom.sizes --blacklist /path/to/custom.blacklist.bed

# Disable BigWig generation for faster processing
nextflow run main.nf --callPeaks --sample_table sample.csv --celltypes celltypes.csv \
  --cistopic.normalize_bigwig false

# Run multiome pipeline with GEX attachment
nextflow run main.nf --callPeaks --inferConsensus --attachGEX --sample_table sample.csv --celltypes celltypes.csv

# Attach GEX data to existing ATAC data
nextflow run main.nf --attachGEX --sample_table updated_sample.csv --celltypes celltypes.csv \
  --atac_adata atac_anndata.csv
```

### Compute Resources
Default resource allocation can be customized in individual config files:
- **CPU cores**: Automatically scaled based on process requirements (1-4 cores)
- **Memory**: Dynamically calculated based on input data size and attempt number
- **Queue**: 'normal' for most processes, 'hugemem' for large datasets
- **Retry strategy**: Up to 5 retries with increased resources on failure

## Key Features
* **Modular design**: Each process is containerized and configurable
* **Scalable**: Memory and CPU requirements auto-adjust based on data size
* **Quality control**: Automatic fragment counting and TSS enrichment analysis
* **Flexible input**: Supports both celltype annotation files and pre-computed pseudobulk peak tables
* **Multiome integration**: Attach GEX data to ATAC data for joint analysis
* **Comprehensive output**: Generates pseudobulks, peaks, consensus peaks, cisTopic objects, AnnData objects, and multiome objects

## Input File Formats

### Sample Table Formats

**Basic sample_table.csv** (for initial peak calling):
```csv
sample_id,path
WS_wEMB13386884,/path/to/cellranger-arc/output/
WS_wEMB13386881,/path/to/cellranger-arc/output/
```

**Updated sample_table.csv** (with fragment counts, generated after peak calling):
```csv
sample_id,path,fragments
WS_wEMB13386884,/path/to/cellranger-arc/output/,730872409
WS_wEMB13386881,/path/to/cellranger-arc/output/,1118846819
```

### Other Input Files

**celltypes.csv** (celltype annotation):
```csv
sample_id,barcode,celltype
WS_wEMB13386884,AGAAGGTGTAATTAGC-1,vasculature
WS_wEMB13386884,GATCGAGCACTTCATC-1,fibroblasts
WS_wEMB13386881,ACAACATGTGATCAGC-1,vasculature
WS_wEMB13386881,GAGCGGTCATGGAGGC-1,fibroblasts
```

**pseudobulk_peaks.csv** (generated after peak calling):
```csv
celltype,fragments,large_peaks,all_peaks,path
vasculature,1000000,5000,5500,/path/to/vasculature_peaks.narrowPeak
fibroblasts,800000,4200,4800,/path/to/fibroblasts_peaks.narrowPeak
```

**atac_anndata.csv** (generated after consensus inference):
```csv
sample_id,path
WS_wEMB13386884,/path/to/WS_wEMB13386884.h5ad
WS_wEMB13386881,/path/to/WS_wEMB13386881.h5ad
```

## Examples of use:
### 1. Perform peak calling
To run peak calling for each celltype you need to specify *sample table* (see `example/sample_table.csv`) and *celltype annotation* (see `example/celltypes.csv`):

```shell
nextflow run main.nf --callPeaks --sample_table ./example/sample_table.csv --celltypes ./example/celltypes.csv
```

This creates `results` directory with the following files:
```
results/
├── log/
│   ├── fragments_celltype_x_sample.csv # fragment counts matrix with shape (n_celltypes, n_samples)
│   ├── fragments_per_celltype.csv # total fragment counts for each celltype
│   ├── pseudobulk.audiovisual_neuroepithelium.log
│   ├── pseudobulk.craniofacial.log
│   └── splitcelltypes.log
├── pseudobulk/
│   ├── fragments/
│   │   ├── audiovisual_neuroepithelium_fragments.tsv.gz
│   │   └── craniofacial_fragments.tsv.gz
│   ├── bigwig/
│   │   ├── audiovisual_neuroepithelium.bw
│   │   └── craniofacial.bw
│   ├── audiovisual_neuroepithelium_peaks.narrowPeak
│   └── craniofacial_peaks.narrowPeak
├── pseudobulk_peaks.csv # contains fragment, peak counts and path to .narrowPeak file for each celltype (see example/pseudobulk_peaks.csv)
└── updated_sample_table.csv # updated sample table which contains fragment counts (see example/updated_sample_table.csv)
```

### 2. Infer consensus peaks and calculate features
To run consensus peak calling and feature calculation you need to specify an **updated sample table** generated on previous step (it is essential to use updated table with fragment counts to set appropriate memory limits for jobs) and **pseudobulk peaks table** generated on previous step with selected celltypes:
```shell
nextflow run main.nf --inferConsensus --sample_table ./results/updated_sample_table.csv --pseudobulk_peaks ./results/pseudobulk_peaks.csv
```

This will create a `consensus_peaks.bed` file, `cisTopic` and `.h5ad` objects for each sample and combined `cisTopic` and `.h5ad` objects for whole dataset:
```
results/
├── consensus_peaks.bed # consensus peaks
├── combined_cistopic_object.pkl
├── combined.h5ad
├── log/
│   └── inferconsensus.log
└── cistopic/
    ├── WS_wEMB13400228/
    │   ├── qc/
    │   ├── WS_wEMB13400228_cistopic_obj.pkl
    │   └── WS_wEMB13400228.h5ad
    └── WS_wEMB13400229/
        ├── qc/
        ├── WS_wEMB13400229_cistopic_obj.pkl
        └── WS_wEMB13400229.h5ad
```

### 3. Attach GEX data to ATAC data
To attach GEX (Gene Expression) data to existing ATAC data for multiome analysis, you need an **updated sample table** (with path pointing to CellRanger-arc output containing both ATAC and GEX data) and **celltype annotation**. You can either provide existing ATAC anndata files or let the pipeline generate them:

Option A - With existing ATAC anndata:
```shell
nextflow run main.nf --attachGEX --sample_table ./example/updated_sample_table.csv --celltypes ./example/celltypes.csv --atac_adata ./results/atac_anndata.csv
```

Option B - Generate ATAC anndata first, then attach GEX:
```shell
nextflow run main.nf --inferConsensus --attachGEX --sample_table ./results/updated_sample_table.csv --celltypes ./example/celltypes.csv --pseudobulk_peaks ./results/pseudobulk_peaks.csv
```

This will create multiome objects (.h5mu) and combined AnnData objects (.h5ad) with both ATAC and GEX data:
```
results/
├── anndata/
│   ├── WS_wEMB13400228/
│   │   ├── WS_wEMB13400228_coupled.h5mu
│   │   └── WS_wEMB13400228_coupled.h5ad
│   └── WS_wEMB13400229/
│       ├── WS_wEMB13400229_coupled.h5mu
│       └── WS_wEMB13400229_coupled.h5ad
└── log/
    └── attachgex.log
```

### 4. Infer consensus peaks and attach GEX in one go
To run consensus peak inference and GEX attachment together:
```shell
nextflow run main.nf --inferConsensus --attachGEX --sample_table ./results/updated_sample_table.csv --celltypes ./example/celltypes.csv --pseudobulk_peaks ./results/pseudobulk_peaks.csv
```

### 5. Perform peak calling, infer consensus peaks and calculate features
To run all steps together you can use the following command:
```shell
nextflow run main.nf --callPeaks --inferConsensus --sample_table ./example/sample_table.csv --celltypes example/celltypes.csv
```

### 6. Complete multiome pipeline (all steps including GEX attachment)
To run the complete pipeline including multiome integration:
```shell
nextflow run main.nf --callPeaks --inferConsensus --attachGEX --sample_table ./example/sample_table.csv --celltypes ./example/celltypes.csv
```