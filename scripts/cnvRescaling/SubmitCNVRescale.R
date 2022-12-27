#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs()
script.dir <- dirname(sub('--file=', '', args[grep('--file=', args)]))
refs.dir <- file.path(script.dir, '..', '..', 'references')

source(file.path(script.dir, 'CNVRescale.R'))

opt_list <- list(make_option('--binned_bed_filename'),
                 make_option('--is_input_control', default = FALSE),
                 # use a flat WCE for cnv_ratios_filename if no WCE exists for APP
                 make_option('--cnv_ratios_filename', default = NULL),
                 make_option('--assembly', default = 'hg19'),
                 # euploid reference filename is used to compare against aneuploid WCEs (in CNAnorm)
                 make_option('--euploid_reference_filename', default = NULL),
                 # if no euploid reference file, then generate one
                 make_option('--simulate_beta', default = 4.81),
                 make_option('--simulate_k', default = 15.7),
                 make_option('--meta_bin_size', default = 50),
                 make_option('--n_windows', default = 10),
                 make_option('--cnv_rescale_output', default = NULL),
                 make_option('--saved_gc_filename', default = file.path(refs.dir, 'hg19_5000_gc.bed')),
                 make_option('--cnv_flag_output_filename', default = NULL),
                 make_option('--cnv_rescale_success_output', default = NULL),
                 make_option('--bypass_cnv_rescaling_step', default = FALSE))

opts <- parse_args(OptionParser(option_list = opt_list))

if(is.null(opts$euploid_reference_filename)){
  euploid_reference_filename <- data.frame('beta' = opts$simulate_beta, 'k' = opts$simulate_k, stringsAsFactors = FALSE)
} else{
  euploid_reference_filename <- opts$euploid_reference_filename
}

is_input_control <- as.logical(opts$is_input_control)

if(is.null(opts$cnv_rescale_output)){
  save_location <- gsub(pattern = '.bed', replacement = '_cnv_rescaled.bed', x = opts$binned_bed_filename)
} else{
  save_location <- opts$cnv_rescale_output
}

# add option to bypass cnv calculation step
if(opts$bypass_cnv_rescaling_step){
  corrected_df <- fread(input = opts$binned_bed_filename, header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)
 } else{
# put this in a try-catch statement
  corrected_df <- tryCatch(expr = RescaleBinnedFileCNV(bin_df_filename = opts$binned_bed_filename, is_input_control = opts$is_input_control,
                                     assembly = opts$assembly, euploid_reference_filename = euploid_reference_filename,
                                     meta_bin_size = opts$meta_bin_size, n_windows = opts$n_windows,
                                     cnv_ratios_filename = opts$cnv_ratios_filename,
                                     saved_gc_filename = opts$saved_gc_filename, save_folder = refs.dir, 
                                     cnv_flag_output_filename = opts$cnv_flag_output_filename), error = function(c) 'error')
}
print(head(corrected_df))
# handle errors outside tryCatch (a little more readable this way)
if(is.data.frame(corrected_df)){
  if(!is.null(opts$cnv_rescale_success_output)){
    write.table(x = 'true', file = opts$cnv_rescale_success_output, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  write.table(x = corrected_df, file = save_location, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
} else{
  if(!is.null(opts$cnv_rescale_success_output)){
    write.table(x = 'false', file = opts$cnv_rescale_success_output, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # if it's an input control, check that cnv_ratios_file, cnv_flag_output_filename exist; if not, write them
  if(opts$is_input_control){
    if(!file.exists(opts$cnv_flag_output_filename)){
      write.table(x = 'error', file = opts$cnv_flag_outpt_filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    if(!file.exists(opts$cnv_ratios_filename)){
      write.table(x = '', file = opts$cnv_ratios_filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
  # write empty cnv_rescaled file
  write.table(x = '', file = save_location, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}





