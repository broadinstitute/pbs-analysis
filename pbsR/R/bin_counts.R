#' Bin genome given chromosome sizes, and count number of reads within each bin
#'
#' This function takes in filepath that specifies chromosome lengths, bins the genome according to bin_size and 
#' gets the count of reads (single-end) or fragments (paired-end) using featureCounts
#' @param bam_file filepath to BAM file. BAM file must be indexed
#' @param chrom_sizes_file filepath to chromosome sizes. Can get from UCSC genome browser. 
#' @param bin_size size of bin to split genome into non-overlapping windows
#' @param paired_end true/false if dataset is paired/single read
#' @param threads number of threads to run featureCounts
#' @export
getBinnedCounts = function(bam_file,
                            chrom_sizes_file,
                            bin_size, 
                            paired_end, 
                            threads = 4){
  
  chrom_sizes = read.table(chrom_sizes_file, header = FALSE, col.names = c("chromosome", "size"))
  
  # Create a GRanges object representing the entire genome
  genome_gr = GenomicRanges::GRanges(
    seqnames = chrom_sizes$chromosome,
    ranges = IRanges::IRanges(start = 1, end = chrom_sizes$size)
  )
  GenomeInfoDb::seqlengths(genome_gr) <- chrom_sizes$size
  
  # Tile the genome according to bin_size 
  tiled_genome = GenomicRanges::tileGenome(seqlengths = GenomeInfoDb::seqlengths(genome_gr), 
                            tilewidth = bin_size, 
                            cut.last.tile.in.chrom = T)
  
  #Create SAF format df for featureCounts
  tiled_genome_df = as.data.frame(tiled_genome)
  tiled_genome_df = tiled_genome_df[,c("GeneID","seqnames","start","end","strand")]
  colnames(tiled_genome_df) = c("GeneID","Chr","Start","End","Strand")
  
  # Use featureCounts to get counts in the tiles
  fc = Rsubread::featureCounts(files = bam_file, 
                     annot.ext = tiled_genome_df,
                     isPairedEnd = paired_end, 
                     nthreads = threads, 
                     isGTFAnnotationFile = F,
                     largestOverlap = T)
  
  return(fc)
}

#TODO: add binnig method for fragment files/ATAC-seq data