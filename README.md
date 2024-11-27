# ATAC pipeline
This pipeline performs peak calling for ATAC data (only cisTopic option is available at the moment)

## Contents of Repo
* `main.nf` - the Nextflow pipeline that runs the whole pipeline
* `modules/` - a collection pipeline's of processes for different components
  * `modules/pycistopic/main.nf` - a collection of processes for cisTopic component
* `workflows/` - a collection pipeline's of workflows for different components
  * `workflows/pycistopic/main.nf` - a collection of processes for cisTopic component
* `bin` - a colection of python scripts
  * `bin/colored_logger.py` - a custom python logger with colored output
  * `bin/make_pseudobulk.py` - a script to make a pseubobulks from fragments file
  * `bin/peak_calling.py` - a script to call peaks for celltypes' pseudobulks
  * `bin/infer_consensus.py` - a script to infer consensus peaks from pseudobulks' peak calling results
  * `bin/create_cistopic.py` - script that creates cisTopic object from consensus `.bed` and fragments files
* `nextflow.config` - the configuration script that controls everything

## Pipeline Arguments
* `--sample_table`: specify a path to .csv file with sample names and path to the CellRanger-arc output dir (see example below)
* `--celltypes`: specify a path .csv file with celltype annotation or a path to pseudobulk_peaks.csv file with selected celltypes for consensus peak calling
* `--callPeaks`: if specified creates pseudobulks for celltypes specified in `--celltypes` for samples in `--sample_table`
* `--inferConsensus`: if specified runs a consensus peak calling and outputs cisTopic object for each sample in `--sample_table`

## Examples of use:
### 1. Perform peak calling
To run peak calling for each celltype you need to specify *sample table* (see `examples/sample_table.csv`) and *celltype annotation* (see `examples/celltype_annotation.csv`):

```shell
nextflow run main.nf --callPeaks --sample_table ./examples/sample_table.csv --celltypes ./examples/celltype_annotation.csv
```

This creates `results` directory with the following files:
```
results/
├── log
│   ├── fragments_celltype_x_sample.csv # fragment counts matrix with shape (n_celltypes, n_samples)
│   ├── fragments_per_celltype.csv # total fragment counts for each celltype
│   ├── pseudobulk.audiovisual_neuroepithelium.log
│   ├── pseudobulk.craniofacial.log
├── narrowPeaks
│   ├── audiovisual_neuroepithelium_peaks.narrowPeak
│   ├── craniofacial_peaks.narrowPeak
├── pseudobulk_peaks.csv # contains fragment, peak counts and path to .narrowPeak file for each celltype (see example/pseudobulk_peaks.csv)
└── updated_sample_table.csv # updated sample table which contains fragment counts (see example/updated_sample_table.csv)
```

### 2. Infer consensus peaks and calculate features
To run consensus peak calling and feature calculation you need to specify an **updated sample table** generated on previous step (it is essential to use updated table with fragment counts to set appropriate memory limits for jobs) and **pseudobulk peaks table** generated on previous step with selected celltypes:
```shell
nextflow run main.nf --inferConsensus --sample_table ./results/updated_sample_table.csv --celltypes ./results/pseudobulk_peaks.csv
```

This will create a `consensus_paeks.bed` file and `cisTopic` objects for each sample:
```
results
├── consensus_peaks.bed # consensus peaks
├── WS_wEMB13400228
│   ├── qc
│   └── WS_wEMB13400228_cistopic_obj.pkl
└── WS_wEMB13400229
    ├── qc
    └── WS_wEMB13400229_cistopic_obj.pkl
```

### 3. Perform peak calling, infer consensus peaks and calculate features
To run all steps together you can use the following command:
```shell
nextflow run main.nf --callPeaks --inferConsensus --sample_table ./examples/sample_table.csv --celltypes examples/celltype_annotation.csv
```