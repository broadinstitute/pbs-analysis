### Calculating the probability of being signal ###
The probability of being signal (PBS) is a value between 0 and 1 assigned to each bin that corresponds to the probability of that bin containing a counts value that corresponds to signal.  Signal can be defined as a peak, such as in H3K27ac or ATAC data, or as a broader, mesa-like region with generally elevated numbers of reads, such as in H3K27me3.  The calculation for PBS can be considered to be a variant of the local false discovery rate (FDR).  False discovery rate, or FDR, is defined as the ratio between the number of false positives and total positives in a dataset; local FDR extends this notion but to a single value.  Thus, for a single value, local FDR will quantify the ratio between the total number of false positives and the total number of positives.  For PBS, we define this ratio based on the shape of the background distribution estimated in the Fitting step of the pipeline, compared to the empirical distribution (i.e. of all of the data).  

Interpreting ChIP-seq data with PBS has several advantages.  First, it considerably simplifies comparing datasets with different read depths, since values are internally normalized compared to the per-dataset background distribution.  Additionally, differentially enriched regions can be easily identified by subtracting two lists of per-bin PBS values.  This increased robustness, however, comes at the cost of lost resolution: PBS is calculated for every 5 kB, and it provides little information about a "strong" vs "very strong" vs "weak" signal, for example.  However, when used in conjunction with other methods, it can provide valuable genomic insight.

### Inputs ###
--binned_bed_filename Binned bed file from binning step of pipeline
--params_df_filename Text file containing parameters of estimated gamma distribution (shape, rate, lambda) from fitting step of pipeline.  
--pbs_filename Output filename for calculated PBS

### Outputs ###
Bed file containing four columns: chr, bin start, bin end, pbs. The genomic positions of the bins matches that of the input binned bed file.  PBS values range from 0 to 1.

Note: additional optional parameters are available to output QC plots to examine the distribution fit to the binned data.
