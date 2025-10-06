# ATAC pipeline
This pipeline performs peak calling for ATAC-seq data using pyCisTopic. It supports pseudobulk creation, peak calling, consensus peak inference, and cisTopic object generation for downstream analysis.

## Contents of Repo
* `main.nf` - the Nextflow pipeline that runs the whole pipeline
* `modules/local/` - a collection of pipeline processes organized by component
  * `modules/local/cistopic/` - cisTopic-related processes (countfragments, splitannotation, pseudobulk, callpeaks, inferconsensus, qualitycontrol, createobject, combineobjects)
  * `modules/local/anndata/` - AnnData-related processes (concat)
* `workflows/` - high-level workflow orchestration
  * `workflows/pycistopic/main.nf` - main pyCisTopic workflow
* `subworkflows/local/` - reusable subworkflow components
  * `subworkflows/local/cistopic_peakcalling/` - peak calling subworkflow
  * `subworkflows/local/cistopic_inferpeaks/` - consensus peak inference subworkflow
* `configs/` - modular configuration files for each process
* `reference/` - reference files (chromosome sizes, blacklists, TSS annotations)
* `nextflow.config` - the main configuration script

## Pipeline Arguments
* `--sample_table`: Path to .csv file with sample names and paths to the CellRanger-arc output directories
* `--celltypes`: Path to .csv file with celltype annotation or path to pseudobulk_peaks.tsv file with selected celltypes for consensus peak calling
* `--callPeaks`: Run peak calling for provided celltypes (needs --sample_table and --celltypes)
* `--inferConsensus`: Run consensus peak calling and create cisTopic objects (needs --sample_table and --celltypes)
* `--output_dir`: Output directory (default: 'results')
* `--help`: Show help message

## Key Features
* **Modular design**: Each process is containerized and configurable
* **Scalable**: Memory and CPU requirements auto-adjust based on data size
* **Quality control**: Automatic fragment counting and TSS enrichment analysis
* **Flexible input**: Supports both celltype annotation files and pre-computed pseudobulk peak tables
* **Comprehensive output**: Generates pseudobulks, peaks, consensus peaks, cisTopic objects, and AnnData objects

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
├── pseudobulk_peaks.tsv # contains fragment, peak counts and path to .narrowPeak file for each celltype (see example/pseudobulk_peaks.csv)
└── updated_sample_table.csv # updated sample table which contains fragment counts (see example/updated_sample_table.csv)
```

### 2. Infer consensus peaks and calculate features
To run consensus peak calling and feature calculation you need to specify an **updated sample table** generated on previous step (it is essential to use updated table with fragment counts to set appropriate memory limits for jobs) and **pseudobulk peaks table** generated on previous step with selected celltypes:
```shell
nextflow run main.nf --inferConsensus --sample_table ./example/updated_sample_table.csv --celltypes ./example/pseudobulk_peaks.tsv
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

### 3. Perform peak calling, infer consensus peaks and calculate features
To run all steps together you can use the following command:
```shell
nextflow run main.nf --callPeaks --inferConsensus --sample_table ./example/sample_table.csv --celltypes example/celltypes.csv
```