#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs()
script.dir <- dirname(sub('--file=', '', args[grep('--file=', args)]))

source(file.path(script.dir, 'ProbabilityBeingSignal.R'))

opt_list <- list(make_option('--binned_bed_filename'),
                 make_option('--params_df_filename'),
                 make_option('--pbs_filename', default = NULL),
                 make_option('--max_dens', default = 10),
                 make_option('--plot_data', default = FALSE),
                 make_option('--title_str', default = NULL),
                 make_option('--return_bin_df', default = FALSE),
                 make_option('--return_plt', default = FALSE),
                 make_option('--plot_distribution_fit_only', default = FALSE))

opts <- parse_args(OptionParser(option_list = opt_list))

if(readLines(con = opts$binned_bed_filename, n = 1) == ''){
  write.table(x = '', file = opts$pbs_filename, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
} else{
  out <- GetProbabilityBeingSignal(bin_df_filename = opts$binned_bed_filename, params_df_filename = opts$params_df_filename,
                                   pbs_filename = opts$pbs_filename, plot_data = opts$plot_data,
                                   title_str = opts$title_str, return_bin_df = opts$return_bin_df, return_plt = opts$return_plt,
                                   plot_distribution_fit_only = opts$plot_distribution_fit_only,
                                   max_dens = opts$max_dens)
}

