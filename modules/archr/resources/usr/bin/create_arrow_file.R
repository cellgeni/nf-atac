#!/usr/bin/env Rscript

####################
# Import libraries #
####################
library(parallel)
library(argparse)
suppressPackageStartupMessages(library(ArchR))

####################
# define functions #
####################


get_parser = function(){

  parser <- ArgumentParser()

  parser$add_argument("--fragments", type = "character",
                      help = "Path to fragment file")

  parser$add_argument("--sample_id", type = "character",
                      help = "Sample id")

  parser$add_argument("--min_counts", type = "integer",
                      help = "minFrags threshold")

  parser$add_argument("--min_tsse", type = "numeric",
                      help = "minTSS threshold")
  
  parser$add_argument("--genome", type = "character",
                      help = "version of genome, only hg38 is available for now")

  parser$add_argument("--nthreads", type = "integer",
                      help = "number of parallel threads")

  parser
}

process_sample = function(fragments,sample_id,min_counts,min_tsse){
  # rhdf5::h5disableFileLocking()
  dir.create(sample_id)
  setwd(sample_id)
  
  ArrowFile = createArrowFiles(
    inputFiles = paste0('../',fragments),
    sampleNames = sample_id,
    outputNames = sample_id,
    minTSS = min_tsse, 
    minFrags = min_counts, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE, 
    force = TRUE
  )

  doubScores <- addDoubletScores(
      input = ArrowFile,
      k = 10, 
      knnMethod = "UMAP", 
      LSIMethod = 1
  )
  
  setwd('..')
}

main = function(){
  parser = get_parser()
  args = parser$parse_args()
  
  jsonlite::write_json(args,'args.json')

  addArchRGenome(args$genome)
  addArchRThreads(args$nthreads)   

  process_sample(fragments = args$fragments,
                 sample_id = args$sample_id,
                 min_counts = args$min_counts,
                 min_tsse = args$min_tsse
  )
  
  data_dir <- args$data_dir
  sample <- args$sample 
  outdir <- args$outdir
}

###############
# Run #########
###############
main()

