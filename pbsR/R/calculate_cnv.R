#' Check for CNVs based on the presence of a bimodal distribution in genome-wide binned counts. 
#' 
#' @param counts_df dataframe with columns chr, start, end, and counts. Generally, the input is the output of pbsR::getMappabilityScore
#' @param meta_bin_size TODO: add DOC
#' @param n_windows TODO: add DOC
testCNVPresence = function(counts_df, meta_bin_size = 50, n_windows = 10, bin_size = 0){
  
  #determine binsize if not user-defined
  if(bin_size == 0){
    bin_size = counts_df$end[1] - counts_df$start[1] + 1
  }
  
  counts_df$rounded_start = floor(counts_df$start/(meta_bin_size*bin_size))*meta_bin_size*bin_size
  rounded_df = counts_df %>% 
    dplyr::group_by(chr, rounded_start) %>% 
    dplyr::filter(counts > 0) %>%
    dplyr::summarise(mean_counts = mean(counts))
  idx_array = rep(x = 1:n_windows, each = ceiling(dim(rounded_df)[1]/n_windows))
  rounded_df$idx = idx_array[1:dim(rounded_df)[1]]
  test_val = rounded_df %>% 
    dplyr::group_by(idx) %>% 
    dplyr::summarise(p_value = dip.test(mean_counts)[[2]],D = dip.test(mean_counts)[[1]])
  
  return(test_val[which.min(test_val$p_value),])
}

#' Given set of genomic coordinates in a dataframe (chr, start, end), get the GC content using BSGenome libraries
#' 
#' @param bed_df dataframe with columns chr, start, end. 
#' @param genome genome name For example, hg19, hg38, mm10 
#' @param bin_size size of genome wide bins. If not user defined, will calculate from bed_df

getGCContent <- function(bed_df, genome, bin_size = 0){
  
  #determine binsize if not user-defined
  if(bin_size == 0){
    bin_size = counts_df$end[1] - counts_df$start[1] + 1
  }
  
  bsgenome = ""
  
  if(genome == "hg19"){
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
      install.packages("BSgenome.Hsapiens.UCSC.hg19")
    bsgenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19  
  }
  else if(genome == "hg38"){
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
      install.packages("BSgenome.Hsapiens.UCSC.hg38")
    bsgenome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 
  }
  else if(genome == "mm10"){
    if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
      install.packages("BSgenome.Mmusculus.UCSC.mm10")
    bsgenome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10 
  }
  else{
    stop("Genome not supported.")
  }
  
  gr = GenomicRanges::makeGRangesFromDataFrame(df = bed_df, seqnames.field = "chr", start.field = "start", end.field = "end")
  seqs <- BSgenome::getSeq(bsgenome, gr)
  
  bed_df$gc = as.numeric(Biostrings::letterFrequency(x = seqs, letters = "GC", as.prob = TRUE)) 
  return(bed_df)
}


#' Calculate CNV ratios
#' 
#' @param counts_df dataframe with columns chr, start, end, and counts. 
#' @param reference_k TODO: add DOC
#' @param reference_beta TODO: add DOC
getCNVWithCNAnorm = function(counts_df, reference_k = 15.7, reference_beta = 4.81){
  
  counts_df = counts_df %>% dplyr::filter(counts > 0)
  
  bin_size = counts_df$end[1] - counts_df$start[1] + 1
  
  reference_df = data.frame('chr' = counts_df$chr, 
                            'start' = counts_df$start, 
                            'end' = counts_df$end,
                            'reference' = stats::rgamma(n = nrow(counts_df), shape = reference_k, rate = reference_beta),
                             stringsAsFactors = FALSE)
  
  covData = dplyr::left_join(x = counts_df, y = reference_df, by = c('chr', 'start', 'end'))
  covData = pbsR:::getGCContent(bed_filename = covData[,c('chr', 'start', 'end')], genome = genome, bin_size = bin_size)
  
  #prep a dataframe for calculating CNV ratios
  df = data.frame("Chr"=covData$chr, "Pos"=covData$start, "Test"=covData$counts, "Norm"=covData$reference, "GC"=covData$gc)
  
  CN = CNAnorm::dataFrame2object(df) %>% 
    CNAnorm::gcNorm(.) %>% 
    CNAnorm::addSmooth(., lambda=7 ) %>% 
    CNAnorm::peakPloidy(., method='closest') %>%
    CNAnorm::validation(.) %>% 
    CNAnorm::addDNACopy(.) %>% 
    CNAnorm::discreteNorm(.) %>% 
    CNAnorm::peakPloidy(., ploidyToTest = 12)
  
  cnv_df = data.frame('chr' = CN@InData@Chr, 
                      'start' = CN@InData@Pos,
                      'end' = CN@InData@Pos + bin_size,
                      'ratio.s.n' = CN@DerivData@ratio.s.n)
  
  cnv_df = cnv_df %>% dplyr::mutate(ratio.s.n = ifelse(test = is.na(ratio.s.n), yes = 1, no = ratio.s.n))
  return(cnv_df)
}

#' Get CNV ratios from a control sample. Control should be binned and processed similarly to treatment sample
#' 
#' @param counts_df dataframe with columns chr, start, end, counts and map_rescaled_counts. 
#' @param meta_bin_size param passed to pbsR::TestCNVPresence
#' @param n_windows param passed to pbsR::TestCNVPresence
getCNVRatioFromControl = function(counts_df, meta_bin_size = 50, n_windows = 10){
  
  bin_size = counts_df$end[1] - counts_df$start[1] + 1
  
  cnv_test_stat = TestCNVPresence(bin_df_filename = counts_df %>% dplyr::filter(chr %in% paste0('chr', 1:22)), 
                                  bin_size = bin_size, 
                                  meta_bin_size = meta_bin_size, 
                                  n_windows = n_windows)
  
  if(cnv_test_stat$p_value < 0.05){
    print('CNVs detected.')
    cnv_df = getCNVWithCNAnorm(counts_df = counts_df)
    # in case a bin doesn't have the same end idx, only join by start
    counts_df <- left_join(x = counts_df, y = cnv_df[,c('chr', 'start', 'ratio.s.n')], by = c('chr', 'start'))
    # correct for NAs
    counts_df$ratio.s.n[is.na(counts_df$ratio.s.n)] = 1
  } 
  else{
    counts_df$ratio.s.n = 1
  }
  return(counts_df)
}

#' Rescale counts by CNV ratios
#' 
#' @param counts_df dataframe with columns chr, start, end, counts and map_rescaled_counts. 
#' @param counts_column which column of the dataframe to apply cnv rescaling on. Default: map_rescaled_counts.
#' @param cnv.ratio vector of ratios to scale counts by. Generally, this vector is computed in pbsR::getCNVRatioFromControl and stored in 'ratio.s.n' column of output 
rescaleCNV = function(counts_df, counts_column = "map_rescaled_counts" , cnv.ratios){
  counts_df$cnv_rescaled_counts <- counts_df[,counts_column]/cnv.ratios
  return(counts_df)
}


