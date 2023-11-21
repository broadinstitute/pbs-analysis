#' Rescale counts in chrX and chrY 
#' 
#' @param counts_df dataframe with columns chr, start, end, and counts. Generally, the input to this function is the output of pbsR::RescaleMappability
#' 
RescaleXY <- function(counts_df){
  if(!('chrY' %in% counts_df$chr)){
    return(counts_df)
  }
  # is this a male sample?
  if('chrX' %in% counts_df$chr & mean(counts_df$map_rescaled_counts[counts_df$chr == 'chrY'])/mean(counts_df$map_rescaled_counts[counts_df$chr == 'chrX']) > 0.5){
    counts_df$map_rescaled_counts[counts_df$chr == 'chrX'] <- 2*counts_df$map_rescaled_counts[counts_df$chr == 'chrX']
    counts_df$map_rescaled_counts[counts_df$chr == 'chrY'] <- 2*counts_df$map_rescaled_counts[counts_df$chr == 'chrY']
  }
  return(counts_df)
}

#' Rescale counts using mappability scores  
#' 
#' @param counts_df dataframe with columns chr, start, end, and counts. Generally, the output of pbsR::getBinnedCounts.
#' @param map_df dataframe with columns chr, start, end, and mappability_score. Look at get("hg19_5000_map_100", asNamespace('pbsR')) to see an example
RescaleMappability <- function(counts_df, map_df, map_threshold = 0.5){
  counts_df = dplyr::left_join(x = counts_df, y = map_df, by = dplyr::join_by('chr', 'start', 'end')) %>%
    dplyr::filter(mappability_score > map_threshold) %>%
    dplyr::mutate(map_rescaled_counts = counts/mappability_score)
  return(counts_df[,c('chr', 'start', 'end', 'counts','map_rescaled_counts')])
}

#' Get mappability score for each bin in binned genome, and scale counts according to mappability score  
#' 
#' @param bam_file filepath to BAM file. BAM file must be indexed
#' @param counts_df  dataframe with columns chr, start, end, and counts. Generally, the input to this function is the output of pbsR::getBinnedCounts
#' @param bin_size size of bin to split genome into non-overlapping windows. If not supplied, will extract from counts_df 
#' @param paired_end true/false if dataset is paired/single read
#' @param map_threshold bins with mappability score lower than this threshold will be filtered out 
#' @export
getMappabilityScore = function(bam_file, counts_df, bin_size = 0, genome, paired_end, map_threshold = 0.5){
  
  #These are the preset available read_lengths for which GEM tracks can be computed [Verify this].   
  if(genome == 'hg38'){
    avail_read_length = c(36, 100)
  }
  else{
    avail_read_length = c(36, 75, 100)
  }
  
  #determine read length in BAM
  read_length = GenomicAlignments::qwidth(GenomicAlignments::readGAlignments(Rsamtools::BamFile(bam_file, yieldSize=1)))
  
  #get closest available read length
  if(paired_end){
    ref_read_length = max(avail_read_length)
  } else{
    ref_read_length = avail_read_length[which.min(abs(avail_read_length - read_length))]
  }
  
  #determine binsize if not user-defined
  if(bin_size == 0){
    bin_size = counts_df$end[1] - counts_df$start[1] + 1
  }
  
  #retrive precomputed mappability bedgraph file. 
  #TODO: add code to handle cases where bedgraph does not exist
  map_file = paste(genome, bin_size, "map",ref_read_length, sep = "_")
  tryCatch( {
    map_df = get(map_file, asNamespace('pbsR'))
  }, error = function(e) {
    print(paste0("Reference map file: ", map_file ," with defined bin size does not exist."))
  })
  
  counts_df = pbsR:::RescaleMappability(counts_df = counts_df, 
                                       map_df = map_df, 
                                       map_threshold = map_threshold)
  counts_df = pbsR:::RescaleXY(counts_df = counts_df)
  return(counts_df)
}
  