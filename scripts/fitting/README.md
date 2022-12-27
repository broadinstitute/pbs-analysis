### Fitting a gamma distribution to background ###
To determine which binned regions contain signal, we first determine which binned regions likely represent background.  To do so, we fit a gamma distribution to the bottom 50th percentile of the data, which we assume exclusively represents background.  This step of the pipeline is implemented by minimizing the Cramer von Mises criterion between an estimated gamma distribution and the truncated distribution of counts using optim, a non-linear optimizer in R.  The output of this step is a file containing the distribution parameters for the estimated gamma distribution.  These include the shape and rate for the background gamma distribution, as well as lambda, which represents the fraction of the data corresponding to background (typically between 0.7 and 1).  Note that this step may take several minutes.

### Inputs ###
binned_bed_filename: name of binned bed file that has been rescaled for mappability and presence of CNVs.
params_output: name of text file containing parameters of estimated gamma distribution

### Outputs ###
Text file containing parameters for estimated gamma distribution
