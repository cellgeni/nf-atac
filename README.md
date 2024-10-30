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
