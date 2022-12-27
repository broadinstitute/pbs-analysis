# binning
This directory contains scripts for binning a BAM file into fixed-width intervals.  

### Inputs:
BAM filename\
Genome assembly\
Bin size

### Output:
BED file with binned reads


## Details
`binning.sh` bins the reads of a BAM file using [IGVtools](https://software.broadinstitute.org/software/igv/igvtools_commandline) count.
It examines the reads in the BAM to determine if it is single or paired end.  This is necessary for read extension.  For single-ended reads, it estimates the actual fragment length to be 200bp. For paired-end reads, it calculates the actual fragment length.  The density of extended reads is computed over windows of the specified width to produce counts per bin.  The resulting Wig file is converted to a BigWig and then to a BED file.

Currently supported genomes: mm9, mm10, hg18, hg19, hg38.

__Data exclusion__:  
IGVtools implicitly excludes any reads with the following SAM flags:  
`0x4` (Unmapped), `0x8` (Mate unmapped), `0x100` (secondary alignment), `0x400` (duplicate), `0x800` (supplementary).  
We also explicitly exclude any reads with MAPQ < 1.  
> Detailed SAM file format information [here](https://samtools.github.io/hts-specs/SAMv1.pdf)  

ENCODE [blacklisted regions](https://www.encodeproject.org/annotations/ENCSR636HFF/) are excluded from the final BED file 

__Flags__:\
`-f` full path for BAM file to be binned. Reads need to be sorted already, or a .bai file needs to exist in the same directory.\
`-w` number of base pairs in each bin [*default=5000*]\
`-g` genome assembly. Options are `mm9` `mm10` `hg18` `hg19` `hg38`. [*default=hg19*]\
`-r` fragment size estimate for single-ended library [*default=200*]\
`-n` output filename prefix. [*default=output*]\
`-s` output filename suffix. [*default=binned*]

Sample input: 
```
./scripts/binning/bin.sh -f /bam_dir/aggregated_aln_006275.bam -n /binned_dir/aggregated_aln_006275_
Output file saved as /binned_dir/aggregated_aln_006275_binned.bed.
```
### Directions for running binning on Broad cluster
*All steps should be run from* ./ProductionCode/ *directory*
1. After cloning repo, run ./tools.sh to download the necessary UCSC executable files to the /tools folder.
2. In binning folder, run ./scripts/binning/bin_on_uger.sh with appropriate flags as listed above.
