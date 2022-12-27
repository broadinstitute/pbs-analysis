#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

RescaleXY <- function(bin_df){
  if(!('chrY' %in% bin_df$chr)){
    return(bin_df)
  }
  # is this a male sample?
  if('chrX' %in% bin_df$chr & mean(bin_df$counts[bin_df$chr == 'chrY'])/mean(bin_df$counts[bin_df$chr == 'chrX']) > 0.5){
    print('Male donor')
    bin_df$counts[bin_df$chr == 'chrX'] <- 2*bin_df$counts[bin_df$chr == 'chrX']
    bin_df$counts[bin_df$chr == 'chrY'] <- 2*bin_df$counts[bin_df$chr == 'chrY']
  }
  return(bin_df)
}

# identify appropriate map track elsewhere
# map_df should be bed file with chr, start, end, score
# bin_df should be bed file with chr, start, end, counts
# this can also be used for rescaling with CNVs (might make sense to use a threshold there too)
RescaleMappability <- function(bin_df, map_df, map_threshold = 0.5){
  bin_df <- bin_df %>% left_join(x = ., y = map_df, by = c('chr', 'start', 'stop')) %>%
    dplyr::filter(score > map_threshold) %>%
    dplyr::mutate(counts = counts/score)
  return(bin_df[,c('chr', 'start', 'stop', 'counts')])
}

# get information about bam filename
# get the read length elsewhere (from a bash script upstream?)
# get paired end elsewhere (from bash script upstream?)
# TODO: should paired end information even be included here?
# map_filename format: ${gen}_${bin}_${read_length}_map.bed
# example map_filename: hg19_5000_map_76.bed
# include data frame in references folder with information about where to download map files
# available read lengths: 36, 75, 100
# for paired end, use read length = 100
GetMapFileParameters <- function(binned_filename, genome = 'hg19', avail_gen = c('mm9', 'mm10', 'hg19', 'hg38'),
                                 map_file_directory = 'references', read_length, paired,
                                 binning_script_filename = file.path('scripts', 'binning', 'bin.sh'),
                                 map_url_filename = file.path('references', 'MapTrackURL.csv'),
                                 avail_read_length = c(36, 75, 100)){
  # which genome assembly is it?
  if(!(genome %in% avail_gen)){
    print(genome)
    stop('genome assembly unavailable', call. = TRUE)
  }
  # no 75 bp available for hg38
  if(genome == 'hg38'){
    avail_read_length <- c(36, 100)
  }
  # what's the closest ref_read_length to the read length of the file?
  if(paired == TRUE){
    ref_read_length <- max(avail_read_length)
  } else{
    ref_read_length <- avail_read_length[which.min(abs(avail_read_length - read_length))]
  }
  # what's the bin size?
  first_line <- fread(input = binned_filename, nrows = 1, data.table = FALSE, stringsAsFactors = FALSE,
                      col.names = c('chr', 'start', 'end', 'counts'))
  bin_size <- first_line$end - first_line$start
  map_filename <- file.path(map_file_directory, paste(genome, bin_size, 'map', ref_read_length, sep = '_')) %>%
                      paste0(., '.bed')
  # create map file if it doesn't exist already
  if(!file.exists(map_filename)){
    print(paste(map_filename, 'does not exist. Creating mappability reference.'))
    # determine appropriate reference genome tiles
    reference_file <- file.path(map_file_directory, paste0(genome, '_', bin_size, '_tiles.bed'))
    map_urls <- fread(input = map_url_filename, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE) %>%
      dplyr::filter(Genome %in% genome)
    # find kmer length that's closest to read_length; if paired, just use longest read length available
    if(paired){
      map_track <- map_urls$url[which.max(map_urls$kmer)] %>% as.character(.)
    } else{
      map_track <- map_urls$url[which.min(abs(map_urls$kmer - read_length))] %>% as.character(.)
    }
    chrom_sizes <- file.path(map_file_directory, paste0(genome, '.chrom.sizes'))
    sys_command <- paste(binning_script_filename,
                         '--mapURL', map_track,
                         '--mapFile', basename(map_track),
                         '--binSize', bin_size,
                         '--referenceFile', reference_file,
                         '--chromSizes', chrom_sizes,
                         '--saveFilename', map_filename)
    system(sys_command)
  }
  return(map_filename)
}

# run this function to bin a mappability file when known that will use this bin size/read length frequently.
# needs to run with a shell script that calls UGER to set appropriate memory, working directory
# NOTE: this function is never run as part of the pipeline
# reference file is hg19_5000_tiles.bed, or hg38_tiles.bed in references folder
# example: GetBinnedMapFile(unbinned_map_filename='~/ProductionCode/tmp.0Cms7p8LJx/k36.umap.bedgraph.gz', reference_filename='~/ProductionCode/references/hg38_5000_tiles.bed', read_length='36', genome='hg38')
GetBinnedMapFile <- function(unbinned_map_filename, reference_filename, genome = 'hg19',
                             bin_size = '5000', read_length = '75', 
                             bin_map_script = '~/ProductionCode/scripts/rescaling/bin_map_file.sh',
                             save_location = '~/ProductionCode/references/'){
  save_filename <- file.path(save_location, paste0(genome, '_', bin_size, '_map_', read_length, '.bed'))
  chrom_sizes <- file.path(save_location, paste0(genome, '.chrom.sizes'))
  sys_command <- paste(bin_map_script,
                       '--mapFile', unbinned_map_filename,
                       '--referenceFile', reference_filename,
                       '--chromSizes', chrom_sizes,
                       '--binSize', bin_size,
                       '--saveFilename', save_filename)
  system(sys_command)
}

