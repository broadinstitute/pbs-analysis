
# Rescaling
Scripts in this folder handle rescaling of binned values to account for mappability of bins.  

### Inputs ###
Binned file\
BAM file\
Genome

### Outputs ###
Rescaled binned bed file\
SNR output file

## Details
### Alignability and mappability ###
At the binning step in the pipeline, reads that cannot be uniquely aligned to the genome are removed.  Alignability of reads is a function of both read length and whether the reads are single-end or paired-end, with longer, paired-end reads being the most alignable and short, single-ended reads the least.  The corresponding value from the genomic perspective is the mappability of a particular base pair, and subsequently of a bin.  Bins with low mappability - i.e. to which should align low alignability reads - will therefore have a an under-reported value of counts-per-bin, as reads with low alignability have been removed.  To avoid underreporting counts-per-bin values for bins with low mappability, we approximated mappability across each bin by averaging the per-bp values reported in GEnome Multitool (Derrien 2012) for read lengths 36 bp or 75 bp, and from Umap for read lengths of 25, 50 or 100 bp (Karimzadeh 2018) across each bin.  Due to the complexity inherent to estimating a per-bp mappability score for paired-end reads, which do not have pre-defined read-lengths, we used mappability values for 100 bp reads for all paired-end data.  Mappability values ranged from 0 (not mappable) to 1 (uniquely mappable).  Bins with mappability values less than 0.5 are excluded, and counts in bins with mappability values between 0.5 and 1 were rescaled by dividing by the mappability score.  Thus, we correct for the fold-reduction in counts-per-bin due to low mappability of reads by scaling up the counts based on the mappability of the bin.

Mappability tracks with 5 kB-binned values for each of the sequence lengths mentioned above (25, 36, 50, 75 and 100) are found in the `./references/` folder.

Additional information regarding the GEnome Multitool (GEM): https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377 
Additional information regarding Umap mappability: https://academic.oup.com/nar/article/46/20/e120/5086676 and https://bismap.hoffmanlab.org/ 

### Handling sex chromosomes ###
In male samples, X and Y chromosomes are each haploid, and therefore the associated values of counts-per-bin will be halved.  To correct for this, we first determine whether the sample is male or female, and then double the counts for X and Y in the case of a male sample.  Note that due to the low mappability of most bins in the Y chromosome, many of the bins will still be ignored.  Additionally, for female samples, some reads are still aligned to the Y chromosome, but these can be considered background noise.

### Calculating a signal to noise ratio ###
The signal to noise ratio (SNR) value calculated here is the ratio between the value of the 99th percentile of counts per bin to the bottom 10th percentile.  Though this value is only a rough approximation, we find that it corresponds fairly well to visual inspection of the data.  We recommend its use alongside the random forest classifier as a quality control metric.

__Flags__\
`--bam_filename` bam file; this is the same as the input from the binnning step, and is needed to calculate the read length\
`--binned_bed_filename` binned file; intermediate output from binning step\
`--genome` genome build; same as from the binning step, needed to determine the appropriate mappability reference\
`--output_filename` string indicating the name of the output file from this step; recommended to use basename used in binning step with different suffix\
`--snr_output_filename` string indicating the name of the output file to store the SNR value; recommended to use basename in binning step with `_snr.txt` suffix



