# CNV rescaling
The scripts in this folder adjust per-bin counts for the presence of copy number variations by calculating a ploidy for each bin in the case of whole cell extract (WCE)/input control, or dividing the bin counts by a previously calculated ratios file.  

### Inputs
Binned bed file from the mappability rescaling step.\
T/F flag indicating whether the file corresponds to an input control.\
Filename where CNV ratios are to be saved if WCE, or where they have been saved if not WCE.\
Option to bypass CNV rescaling step (e.g. if the datatype is ATAC, or some other that does not have an input control).

### Outputs ###
Binned bed file rescaled for presence of CNVs, stored at cnv_rescale_output.\
If the binned bed file is an input control, a binned bed file with a value at each bin indicating the ploidy, and which will be used in subsequently rescaling epitope files.\
Text file indicating whether CNVs were detected in the binned file, stored at cnv_flag_output_filename.\
Text file indicating whether this step of the pipeline ran correctly, stored in cnv_rescale_success_output.

## Details
### Copy number variations ###
Non-diploid chromosomal regions, most often found in tumor samples or cancer cell lines, can lead to artificial increases or decreases in bin counts.  Copy number variations, or CNVs, are most easily detected in input control datasets, as ploidy differences can be obscured by peaks or true signal, particularly in the case of broad histone marks.  Therefore, correcting for CNVs most rigorously requires a matching input control for any epitope track; due to tumor heterogeneity, the input control should match the epitope biosample as closely as possible.  CNV detection and ploidy estimation steps are calculated on input control tracks; the ratios resulting from the ploidy estimation step are used to correct the epitope tracks.

### Pipeline components
We first determine whether a sample contains CNVs by splitting the entire genome of a WCE into 20 parts, and then using Hartigans’ dip test to determine multimodality of any of the subsets.  By splitting the genome into a number that is not a multiple of the number of chromosomes, we increase our ability to detect CNVs that are due to loss of entire chromosomes; the large number of divisions improves our sensitivity to losses of relatively small regions (Figure?).  Additionally, this CNV detection step is very rapid (< 1 second) and saves a considerable amount of time by not unnecessarily calculating the ratios required to correct for aneuploidy.   If CNVs are detected, we use the CNAnorm package (Gusnanto, 2012) to correct for them.  Briefly, CNAnorm fits a multimodal normal distribution to a histogram of counts, defines diploid as the mode with the highest frequency of counts, and outputs a value for each bin indicating its fold-change from diploid (e.g. a fold change of 1.5 indicates triploid, and 0.5 indicates haploid).  These fold-change values are saved as a bed file, and are used to rescale the input control file, as well as any epitope tracks associated with the biosample of the input control.  

### Bypassing the pipeline ###
The CNV pipeline relies on the presence of an input control or whole cell extract (WCE).  If none is available, a bed file with the same bins locations as the current file being analyzed but with a fourth column with values equal to one can be used as a dummy.  Alternatively, the flag `--bypass_cnv_rescaling_step` can be set to TRUE, and then the final binned file will be the same as the file output from the mappability rescaling step.  Both of these options yield the same result, and implementing either of them is recommended only in cases when one is certain that the dataset does not contain CNVs.    

__Flags__\
`--binned_bed_filename` file output from rescaling step of the pipeline (the step immediately prior to this one)\
`--is_input_control` whether the file is an input control\
`--cnv_ratios_filename` location of the CNV ratios bed file\
`--assembly` genome assembly\
`--cnv_rescale_output` filename where file corrected for CNVs is to be stored if WCE, or where it is stored if not WCE\
`--saved_gc_filename` file with pre-calculated values for GC content.  If this file is not provided, it will be calculated based on a reference and assembly\
`--cnv_flag_output_filename` filename for where information regarding whether CNVs have been detected in the binned file\
`--cnv_rescale_success_output` filename indicating where information is stored regarding whether the CNV step of the pipeline completed successfully\
`--bypass_cnv_rescaling_step` if TRUE, allows to bypass this step entirely, and the final binned file will be the same as the one output by the mappability rescaling step.

