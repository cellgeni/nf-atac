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

  parser$add_argument("--arrow", type = "character",nargs='+',
                      help = "Path to arrow file")

  parser$add_argument("--out_base", type = "character",
                      help = "base name for output (mtx cpm matrix, obs and var csv)")
  
  parser$add_argument("--genome", type = "character",
                      help = "version of genome, only hg38 is available for now")

  parser$add_argument("--nthreads", type = "integer",
                      help = "number of parallel threads")

  parser
}


get_gene_scores = function(arrows){
    prj = ArchRProject(
        ArrowFiles = arrows, 
        outputDirectory = 'tmp',
        copyArrows = FALSE
    )

    gene_scores = getMatrixFromProject(prj,useMatrix = 'GeneScoreMatrix')
    colnames(gene_scores) = sub('#',':',colnames(gene_scores))
    gene_scores
}

main = function(){
  parser = get_parser()
  args = parser$parse_args()

  jsonlite::write_json(args,'args.json')
  
  addArchRGenome(args$genome)
  addArchRThreads(args$nthreads)

  gene_scores = get_gene_scores(args$arrow)
  
  var = rowData(gene_scores)
  rownames(var) = paste0(var$seqnames,':',var$start,'-',var$end)
  write.csv(gene_scores@colData,paste0(args$out_base,'_obs.csv'))
  write.csv(var,paste0(args$out_base,'_var.csv'))
  Matrix::writeMM(t(gene_scores@assays@data$GeneScoreMatrix),paste0(args$out_base,'.mtx'))
}

###############
# Run #########
###############
main()