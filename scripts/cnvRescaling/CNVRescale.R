#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(diptest))
# packages for CNANorm
suppressPackageStartupMessages(library(CNAnorm))
#suppressPackageStartupMessages(library(Rsubread))
#suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
# this package does not exist in the RStudio on the Broad server; moved load package to RescaledBinnedFileCNV function
#suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
# this package does not exist in all moved load package to RescaledBinnedFileCNV function
# suppressPackageStartupMessages(library(biovizBase))
# biovizBase is for calculating GC content

# use this function to rescale binned files for CNVs; maybe rethink logic based on LIMS queries
#
RescaleBinnedFileCNV <- function(bin_df_filename, is_input_control = FALSE, cnv_ratios_filename,
                                 assembly = 'hg19', euploid_reference_filename,
                                 meta_bin_size = 50, n_windows = 10, saved_gc_filename = './references/hg19_5000_gc.bed',
                                 save_folder = 'references', cnv_flag_output_filename = NULL,
                                 cnv_rescale_success_filename = NULL) {
  bin_df <- fread(input = bin_df_filename, col.names = c('chr', 'start', 'end', 'counts'),
                  stringsAsFactors = FALSE, data.table = FALSE)
  bin_size <- bin_df$end[1] - bin_df$start[1]
  if(is_input_control){
    # only test somatic chromosomes (leave out X and Y)
    cnv_flag <- GetCNVFlagIdx(bin_df_filename = bin_df %>% dplyr::filter(chr %in% paste0('chr', 1:22)), bin_size = bin_size, meta_bin_size = meta_bin_size,
                              n_windows = n_windows)
    if(!is.null(cnv_flag_output_filename)){
      write(x = ifelse(test = cnv_flag$p_value < 0.05, yes = "true", no = "false"), file = cnv_flag_output_filename)
    }
    if(cnv_flag$p_value < 0.05){
      print('CNVs detected.')
      out <- GetCNVWithCNAnorm(bin_size = bin_size, sample_filename = bin_df_filename,
                               reference_filename = euploid_reference_filename,
                               save_ratios_filename = cnv_ratios_filename, assembly = assembly, return_cnv_df = TRUE,
                               saved_gc_filename = saved_gc_filename, save_folder = save_folder)
      # in case a bin doesn't have the same end idx, only join by start
      bin_df <- left_join(x = bin_df, y = out[,c('chr', 'start', 'ratio.s.n')], by = c('chr', 'start'))
      # correct for NAs
      bin_df$ratio.s.n[is.na(bin_df$ratio.s.n)] <- 1
      bin_df$counts <- bin_df$counts/bin_df$ratio.s.n
      return(bin_df[,c('chr', 'start', 'end', 'counts')])
    } else{
      cnv_df <- bin_df[,c('chr','start', 'end')]
      cnv_df$ratio <- 1
      write.table(x = cnv_df, file = cnv_ratios_filename, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
      return(bin_df)
    }
  }
  cnv_ratios <- fread(input = cnv_ratios_filename, col.names = c('chr', 'start', 'end', 'ratio'),
                      stringsAsFactors = FALSE, data.table = FALSE)
  if(!is.null(cnv_flag_output_filename)){
    write(x = ifelse(test = mean(cnv_ratios$ratio == 1) < 1, yes = "true", no = "false"), file = cnv_flag_output_filename)
  }
  bin_df <- inner_join(x = bin_df, y = cnv_ratios[, c('chr', 'start', 'ratio')], by = c('chr', 'start'))
  # correct for NAs
  bin_df$ratio[is.na(bin_df$ratio)] <- 1
  bin_df$counts <- bin_df$counts/bin_df$ratio
  return(bin_df[,c('chr', 'start', 'end', 'counts')])
}

# this function checks bed file for CNVs based on the presence of a bimodal distribution
# input: bin_df_filename (either a filename or a data frame) that has been binned into 5 kB bins & rescaled for mappability
# expect standard bed format (4 cols: chr, start, stop, counts)
# meta_bin_size: how many bins should be averaged together
# n_windows: how many windows to split genome to detect bimodality?
GetCNVFlagIdx <- function(bin_df_filename, meta_bin_size = 50, n_windows = 10, plot_data = FALSE,
                          bin_size = 5000, return_rounded_df = FALSE, save_flag_location = NULL){
  if(!is.data.frame(bin_df_filename)){
    bin_df <- fread(input = bin_df_filename, col.names = c('chr', 'start', 'end', 'counts'),
                         stringsAsFactors = FALSE, data.table = FALSE)
  } else{
    bin_df <- bin_df_filename
  }
  bin_df$rounded_start <- floor(bin_df$start/(meta_bin_size*bin_size))*meta_bin_size*bin_size
  bin_df$counts[bin_df$counts == 0] <- NA
  rounded_df <- bin_df %>% dplyr::group_by(chr, rounded_start) %>%
    dplyr::summarise(mean_counts = mean(counts, na.rm = TRUE))
  idx_array <- rep(x = 1:n_windows, each = ceiling(dim(rounded_df)[1]/n_windows))
  rounded_df$idx <- idx_array[1:dim(rounded_df)[1]]
  test_val <- rounded_df %>% dplyr::group_by(idx) %>% dplyr::summarise(p_value = dip.test(mean_counts)[[2]],
                                                                       D = dip.test(mean_counts)[[1]])
  if(plot_data){
    plt <- ggplot(rounded_df, aes(x = mean_counts, y = ..density..)) + geom_histogram(binwidth = 0.05) + facet_wrap(~idx) +
      theme_bw(base_size = 16)
    print(plt)
  }
  if(!is.null(save_flag_location)){
    write.table(x = test_val, file = save_flag_location, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  }
  if(return_rounded_df){
    return(rounded_df)
  } else{
    return(test_val[which.min(test_val$p_value),])
  }
}

# need a reference with known normal ploidy
# use same reference always?
GetCNVWithCNAnorm <- function(bin_size = 5000, saved_gc_filename = NULL, assembly = 'hg19',
                              sample_name = 'sample', sample_filename,
                              reference_name = 'reference', reference_filename,
                              save_ratios_filename = NULL, save_cnv_object_filename = NULL,
                              save_folder = 'references', return_cnv_df = FALSE){
  sample_df <- fread(input = sample_filename, col.names = c('chr', 'start', 'end', sample_name), stringsAsFactors = FALSE,
                     data.table = FALSE)
  # remove bins with zeros (this is mainly to avoid issues with female samples and noisy chromosome Y)
  sample_df <- sample_df[sample_df[,sample_name] > 0,]
  # ok for reference_df to have more rows than sample because joining with left_join
  if(is.data.frame(reference_filename)){
    reference_df <- data.frame('chr' = sample_df$chr, 'start' = sample_df$start, 'end' = sample_df$end,
                               'V1' = rgamma(n = dim(sample_df)[1], shape = reference_filename$k, rate = reference_filename$beta),
                               stringsAsFactors = FALSE)
    names(reference_df)[4] <- reference_name
  } else{
    reference_df <- fread(input = reference_filename, col.names = c('chr', 'start', 'end', reference_name),
                        stringsAsFactors = FALSE, data.table = FALSE)
  }
  covData <- left_join(x = sample_df, y = reference_df, by = c('chr', 'start', 'end'))
  # correct for GC content in bins; GCcontent is in biovizBase; may need to be reloaded?
  if(file.exists(saved_gc_filename)){
    gc_df <- fread(input = saved_gc_filename, col.names = c('chr', 'start', 'end', 'gc'))
  } else{
    print('Calculating GC content...')
    gc_df <- GetGCContent(bed_filename = covData[,c('chr', 'start', 'end')], save_folder = save_folder, assembly = assembly, return_gc = TRUE)
  }
  covData <- left_join(x = covData, y = gc_df, by = c('chr', 'start', 'end'))
  if("dplyr" %in% (.packages())){
    detach("package:dplyr", unload=TRUE)
  }
  suppressPackageStartupMessages(library(magrittr))
  df <- data.frame(Chr=covData$chr, Pos=covData$start, Test=covData[,sample_name], Norm=covData[,reference_name],
                   GC=covData$gc)
  # this is the meat of the function
  print('Calculating CNV ratios...')
  CN <- dataFrame2object(df) %>% gcNorm(.) %>% addSmooth(., lambda=7 ) %>% peakPloidy(., method='closest') %>%
          validation(.) %>% addDNACopy(.) %>% discreteNorm(.) %>% peakPloidy(., ploidyToTest = 12)
 # reattach dplyr
  suppressPackageStartupMessages(library(dplyr))
  if(!is.null(save_cnv_object_filename)){
    saveRDS(object = CN, file = save_cnv_object_filename)
#     function_call <- capture.output(print(match.call()))
#     metadata_str <- paste(function_call, '\nsample_name =', sample_name, '\nsample_filename = ', sample_filename,
#                           '\nreference_name =', reference_name, '\nreference_filename =', reference_filename,
#                           '\nsave_filename =', save_filename)
#     metadata_filename <- gsub(pattern = '.rds', replacement = '_fun_params.txt', x = save_filename)
#     write(x = metadata_str, file = metadata_filename)
  }
  if(!is.null(save_ratios_filename)){
    cnv_df <- data.frame('chr' = CN@InData@Chr, 'start' = CN@InData@Pos,
                         'end' = CN@InData@Pos + bin_size,
                         'ratio.s.n' = CN@DerivData@ratio.s.n)
    cnv_df <- cnv_df %>% dplyr::mutate(ratio.s.n = ifelse(test = is.na(ratio.s.n), yes = 1, no = ratio.s.n))
    write.table(x = cnv_df, file = save_ratios_filename, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  if(return_cnv_df){
    return(cnv_df)
  }
}

# use this function to make and save a GC content reference file
# pass in a bed file to get genomic locations, and then save in references folder; bed file has to be >= 3 cols
# saved file format: gc_content_<bin_size>_<assembly>
GetGCContent <- function(bed_filename, save_folder = 'references', assembly = 'hg19', bin_size = 5000,
                         return_gc = FALSE){
  if("dplyr" %in% (.packages())){
    detach("package:dplyr", unload=TRUE)
  }
  library(magrittr)
  if(is.data.frame(bed_filename)){
    binned_df <- bed_filename
  } else{
    binned_df <- fread(input = bed_filename, select = 1:3, stringsAsFactors = FALSE,
                       data.table = FALSE, col.names = c('chr', 'start', 'end'))
  }
  suppressPackageStartupMessages(library(biovizBase))
  # add ones to each start
  binned_df$start <- binned_df$start + 1
  if(grepl(pattern = 'hg19', x = assembly)){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    gcContent <- GCcontent(BSgenome.Hsapiens.UCSC.hg19, as(binned_df, 'GRanges'))[,1]
  } else if(grepl(pattern = 'hg38', x = assembly)){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    gcContent <- GCcontent(BSgenome.Hsapiens.UCSC.hg38, as(binned_df, 'GRanges'))[,1]
  } else{
    stop('genome assembly unavailable', call. = TRUE)
  }
  binned_df$start <- binned_df$start - 1
  binned_df$gc <- gcContent
  if(!is.null(save_folder)){
    save_filename <- file.path(save_folder, paste0(assembly, '_', bin_size, '_gc.bed'))
    write.table(x = binned_df, file = save_filename, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  library(dplyr)
  if(return_gc){
    return(binned_df)
  }
}

# use this function to generate a euploid reference file based on fit parameters of a whole cell extract track
MakeEuploidReference <- function(genome_tiles_filename, params_df_filename, save_filename){
  bin_df <- fread(input = genome_tiles_filename, col.names = c('chr', 'start', 'end'), select = 1:3,
                  stringsAsFactors = FALSE, data.table = FALSE)
  params_df <- fread(input = params_df_filename, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
  sim_counts <- rgamma(n = 1:nrow(bin_df), shape = params_df$k, rate = params_df$beta)
  bin_df$counts <- sim_counts
  write.table(x = bin_df, file = save_filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
}


