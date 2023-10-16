## To install. Might also need to import margrittr and forEach
```remotes::install_github("broadinstitute/pbs-analysis/pbsR@sw-pbs")```

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

