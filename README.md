# pbs-analysis
This repository contains a pipelne for binning, estimating the background distribution, and calculating the probability of being signal for ChIP-seq and related datatypes.  

### Overall inputs:
- BAM file
- Genome assembly name
- Bin size (minimum 5 kB)

### Overall outputs:
- BED file with binned reads that are adjusted for mappability and copy number variations
- Parameters file with estimated values defining a gamma distribution describing the background of the bed file.  This file can be used to calculate per-bin p-values relative to the background distribution.
- BED file with probability of being signal (PBS) values for each bin.

### Setup
Clone repo:
```
git clone git@github.com:broadinstitute/pbs-analysis.git
```

Install dependencies:
```
cd pbs-analysis
Rscript ./scripts/install.R
conda install -c bioconda samtools=1.11 ucsc-wigtobigwig=357 ucsc-bigwigaverageoverbed=366 igvtools=2.14.1
```

### Process a single file
For a non-input control file that does not have CNVs (simplest case):
```
./scripts/run_pipeline.sh -f <bam_filename> --is_input_control FALSE --params_output params.txt \
  --pbs_output pbs.bed --cnv_rescale_success_output cnv_rescale.txt --bypass_cnv_rescaling_step TRUE
```

For a file that has suspected CNVs:
1. Process input control file to create the `cnv_ratios.txt`:
```
./scripts/run_pipeline.sh -f <bam_filename> --is_input_control TRUE --params_output params_input_control.txt \
  --pbs_output pbs_input_control.bed --cnv_rescale_success_output cnv_rescale_input_control.txt \
  --bypass_cnv_rescaling_step FALSE --cnv_ratios_filename cnv_ratios.txt
```

>Note that if the fourth column in `cnv_ratios.txt` is all `1`, then no CNVs were detected in the input control.

2. Process an epitope file, using `cnv_ratios.txt` as the file for rescaling:
```
./scripts/run_pipeline.sh -f <bam_filename> --is_input_control FALSE --params_output params.txt \
  --pbs_output pbs.bed --cnv_rescale_success_output cnv_rescale.txt --bypass_cnv_rescaling_step FALSE \
  --cnv_ratios_filename cnv_ratios.txt
```


## Producing binned files
`./scripts/binning/binning.sh` bins the reads of a BAM file using IGVtools count (https://software.broadinstitute.org/software/igv/igvtools_commandline). 
Determines whether a BAM is single or paired-end.  Calculates the actual fragment length for paired reads, and estimates 200 bp for single-ended reads.  Reads must have minMapQuality > 1 to be included.\
Currently supported genomes: mm9, mm10, hg18, hg19, hg38.\
ENCODE blacklisted regions are excluded from the final BED file (https://www.encodeproject.org/annotations/ENCSR636HFF/)

Flags:\
`-f` full path for BAM file to be binned. Reads need to be sorted already, or a .bai file needs to exist in the same directory.\
`-w` number of base pairs in each bin [default=`5000`]\
`-g` genome assembly. Options are mm9 mm10 hg18 hg19 hg38. [default=`hg19`]\
`-n` output filename prefix. [default=`output`]\
`-s` output filename suffix. [default=`binned`]

## Rescaling binned files
`./scripts/rescaling/SubmitRescaleBinnedFiles.R` adjusts values in each bin to account for the fact that IGVtools count excludes reads with minMapQuality < 1 using mappability scores for each 5 kB region based on the alignability of that region.  Alignability is a function of the uniqueness of the sequence, as well as the genome build and the read length (longer reads are more alignable).  Bins with lower alignability will have fewer uniquely aligned reads assigned to them; the rescaling step adjusts the counts accordingly.\
Mappability scores are well-defined for single-ended reads, but are computationally expensive to calculate from scratch; to save computation time, read lengths are rounded to the nearest available value (36, 50, 75, 100).  The mappability of paired-end reads is approximated as that of a 100-bp single-ended read for the correct genomic build.  More information can be found here: https://bismap.hoffmanlab.org/ and https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377.
Counts in bins with mappability scores > 0.5 are divided by the mappability score.  Bins with mappability scores < 0.5 are removed from the bed file.

Flags:\
`--bam_filename` full path for BAM file to be binned. Reads need to be sorted already, or a .bai file needs to exist in the same directory.\
`--binned_bed_filename` filename output from binning step\
`--genome` genome assembly. Options are mm9 mm10 hg18 hg19 hg38. [default=`hg19`]\
`--save_filename` where to save file. [default=`binned_bed_basename_map_rescaled.bed`]

## Correcting bin counts for copy number variations (CNVs)
`./scripts/cnvRescaling/SubmitCNVRescale.R` outputs a bed file of bins that have been rescaled based on copy number variations, discussed below.  If the file is an input control, then the script will also produce a file named by `cnv_ratios_filename`, which can be used to rescale epitope tracks with the same biosample.\
A CNV, or copy number variation, describes a genetic event in which the number of copies of a region is either greater or less than two.  CNVs must be considered, for example, to avoid falsely identifying a region as enriched for a histone mark when the reason for the greater signal is that the sample has twice the amount of chromatin due to a chromosomal arm being replicated.\
CNVs can only be accurately detected in input control (whole cell extract) tracks originating from the same biosample as an epitope track (same tissue type is not sufficient).  Therefore, certain technologies such as ATAC-seq, are not suitable for this step.

Flags:\
`--binned_bed_filename` filename output from binning step\
`--is_input_control` is the file an input control/whole cell extract? [default=`FALSE`]\
`--cnv_ratios_filename` filename with the ploidy ratios; this is a required parameter if `is_input_control = FALSE`, and should be matched based on biosample [default=`NULL`]\
`--assembly` [default=`hg19`]\
`--euploid_reference_filename` required for calculating the CNV ratios for an input control; can be any euploid track with matching bin size and genome assembly to the input control being analyzed; for 5 kB bins and hg19, a file exists in the references folder and a filename need not be provided [default=`NULL`]\
`--meta_bin_size` [default=`50`]\
`--n_windows default` [default=`10`]\
`--save_cnv_rescale_filename` [default=`NULL`]

## Estimating the background distribution
`./scripts/fitting/SubmitFitDistributionWithCVM.R` outputs the parameters of the background gamma distribution of bin counts.  The output `params.txt` file is used in calculating PBS in the subsequent steps.  

Flags:\
`--binned_bed_filename` filename output from binning step\
`--save_params_filename` where to save parameters; if no location given, will be saved as base_filename_params.txt, where base_filename is the binned_bed_filename without the .bed.

## Calculating PBS
`./scripts/pbs/SubmitProbabilityOfBeingSignal.R` outputs a bed file with PBS values calculated at each bin.

Flags:\
`--binned_bed_filename` Binned bed file from binning step of pipeline.\
`--params_df_filename` Text file containing parameters of estimated gamma distribution (shape, rate, lambda) from background estimation step of pipeline.\
`--pbs_filename` Output filename for calculated PBS

