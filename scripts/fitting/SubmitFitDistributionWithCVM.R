#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs()
script.dir <- dirname(sub('--file=', '', args[grep('--file=', args)]))

source(file.path(script.dir, 'FitDistributionWithCVM.R'))

opt_list <- list(make_option('--binned_bed_filename'),
                 make_option('--params_output'),
                 make_option('--max_dens', default = 10),
                 make_option('--plot_data', default = FALSE),
                 make_option('--plot_terra', default = FALSE),
                 ### these next parameters don't often change
                 make_option('--lambda_range', default = seq(0.6, 0.9, by = 0.1)),
                 make_option('--length_out', default = 3),
                 make_option('--fix_weight', default = TRUE),
                 make_option('--weight_value', default = 500),
                 make_option('--bin_width', default = 0.05))

opts <- parse_args(OptionParser(option_list = opt_list))
binned_bed_filename <- as.character(opts$binned_bed_filename)
save_filename <- opts$params_output
# read first line of bin_df; if the file is empty, write empty files for params_df
if(readLines(con = binned_bed_filename, n = 1) == ''){
  write.table(x = '', file = save_filename, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
} else{
  bin_df <- fread(input = binned_bed_filename, col.names = c('chr', 'start', 'end', 'counts'), data.table = FALSE,
                  stringsAsFactors = FALSE)
  
  params_df <- GetDistributionParametersWithOptim(working_df = bin_df, lambda_range = opts$lambda_range,
                                                  length_out = opts$length_out, fix_weight = opts$fix_weight,
                                                  weight_value = opts$weight_value, plot_data = opts$plot_data,
                                                  plot_terra = opts$plot_terra, bin_width = opts$bin_width, max_dens = opts$max_dens)  
  write.table(x = params_df, file = save_filename, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}


