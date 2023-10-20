## To install and use
```
remotes::install_github("broadinstitute/pbs-analysis/pbsR@sw-pbs")
library(pbsR)

#known bug:  might also need to import margrittr and forEach
library(foreach)
library(margrittr)
```

## To run the package end-to-end:
```
#pass fit_regions_df if you want to do compartment/region correction
pbs_obj = PBS(bin_df, theta = 0.5, fit_regions_df = NULL)
```

## If you want to run step-by-step:
```
#get distribution parameters
params_df = pbsR:::getDistributionParametersWithOptim(working_df = bin_df, theta = theta)

#get PBS scores
pbs_obj = pbsR:::getProbabilityBeingSignal(bin_df = bin_df, params_df = params_df, theta = theta)
```

## To extract results or produce plots
```
#get PBS plot
pbsR::getPBSPlot(pbs_obj)

#write PBS results to file
pbsR::writePBS(pbs_obj, "path_to_output.tsv")

```
