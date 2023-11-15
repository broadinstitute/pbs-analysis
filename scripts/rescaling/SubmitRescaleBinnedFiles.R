#!/usr/bin/env Rscript

options(scipen = 999)

library(optparse)

# infer script directory from command arguments
args <- commandArgs()
script.dir <- dirname(sub('--file=', '', args[grep('--file=', args)]))
refs.dir <- file.path(script.dir, '..', '..', 'references')
bin_map_script <- file.path(script.dir, 'bin_map_file.sh')

source(file.path(script.dir, 'RescaleBinnedFiles.R'))

# parse input arguments
opt_list <- list(make_option('--bam_filename'),
                 make_option('--binned_bed_filename'),
                 make_option('--genome', default = 'hg19'),
                 make_option('--output_filename', default = NULL),
                 make_option('--snr_output_filename', default = NULL))

opts <- parse_args(OptionParser(option_list = opt_list))

# first determine whether files are paired or single end and read length
sys_command <- paste(file.path(script.dir, 'read_stats.sh'),  '-f', opts$bam_filename)
read_stats <- system(sys_command, intern = TRUE)

paired <- read_stats[1] == 1
#paired <- read_stats[1] == 2
read_length <- as.numeric(read_stats[2])

# identify appropriate mappability file
map_filename <- GetMapFileParameters(binned_filename = opts$binned_bed_filename, genome = opts$genome,
                                     read_length = read_length, paired = paired,
                                     map_url_filename = file.path(refs.dir, 'MapTrackURL.csv'),
                                     map_file_directory = refs.dir,
                                     binning_script_filename = bin_map_script)

bin_df <- fread(input = opts$binned_bed_filename, col.names = c('chr', 'start', 'stop', 'counts'),
                data.table = FALSE, stringsAsFactors = FALSE)

# rescale by mappability
map_df <- fread(input = map_filename, col.names = c('chr', 'start', 'stop', 'score'),
                data.table = FALSE, stringsAsFactors = FALSE)
bin_df <- RescaleMappability(bin_df = bin_df, map_df = map_df)

# if X and Y chromosomes haven't been removed with mappability file, rescale them.
if('chrX' %in% bin_df$chr){
  bin_df <- RescaleXY(bin_df = bin_df)
}

# save file
if(is.null(opts$output_filename)){
  output_filename <- gsub(pattern = 'binned.bed', replacement = 'map_scaled.bed', x = opts$binned_bed_filename)
} else{
  output_filename <- opts$output_filename
}

write.table(x = bin_df, file = output_filename, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

# get signal to noise ratio filename
if(is.null(opts$snr_output_filename)) {
  snr_output_filename <- gsub(pattern = 'binned.bed', replacement = 'snr.txt', x = opts$binned_bed_filename)
} else {
  snr_output_filename <- opts$snr_output_filename
}

# write signal to noise ratio; -1 if denominator is 0
quantiles <- quantile(x = bin_df$counts, probs = c(0.1, 0.99))
x <- if (quantiles[1] == 0) -1 else quantiles[2] / quantiles[1]
write(x = x, file = snr_output_filename)
