#!/usr/bin/env Rscript

# Generate reference files of bins
# Expects genome assembly first, bin width second
library('rtracklayer')

args <- commandArgs()

script.dir <- dirname(sub('--file=', '', args[grep('--file=', args)]))
refs.dir <- file.path(script.dir, '..', '..', 'references')

trailArgs <- commandArgs(trailingOnly=TRUE)

genomeInfo.file <- sprintf(file.path(refs.dir, '%s.rds'), trailArgs[1])
blacklist.file <- sprintf(file.path(refs.dir, '%s_blacklist.bed'), trailArgs[1])

genomeInfo <- tryCatch(
  {
    readRDS(genomeInfo.file)
  }, warning=function(cond){
    message("ERROR: Genome assembly not available")
    stop()
  })

blacklist <- tryCatch(
  {
    import.bed(blacklist.file)
  }, warning=function(cond){
    message("Warning: No blacklist available for this assembly")
    return(NULL)
  })

# Determine species
species <- substr(trailArgs[1], 1, 2)
if (species == 'hg'){
  contigs <- paste0('chr', c(seq(22), 'X', 'Y'))
} else if (species == 'mm') {
  contigs <- paste0('chr', c(seq(19), 'X', 'Y'))
}

bins <- tileGenome(genomeInfo[contigs], tilewidth=as.integer(trailArgs[2]), cut.last.tile.in.chrom=T)
bins$name <- seq(length(bins))

# Exclude blacklist
if (!is.null(blacklist)){
  bins <- bins[bins %outside% blacklist]
}

export.bed(bins, sprintf(file.path(refs.dir, '%s_%s_tiles.bed'), trailArgs[1], trailArgs[2]))
