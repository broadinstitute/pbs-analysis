#' @export
PBS <-
  function(bin_df,
           theta = 0.5,
           plot_pbs = TRUE,
           fit_regions_df = NULL) {
    
    if (!is.data.frame(bin_df)) {
      bin_df <-
        read.table(
          file = bin_df,
          stringsAsFactors = FALSE,
          col.names = c("chr", "start", "end", "counts")
        )
    } else {
      colnames(bin_df) <- c("chr", "start", "end", "counts")
    }
    
    params_df <-
      pbsR:::getDistributionParametersWithOptim(working_df = bin_df, theta = theta)
    pbs_obj = pbsR:::getProbabilityBeingSignal(bin_df = bin_df, params_df = params_df, theta = theta)
    
    if (!is.null(fit_regions_df)) {
      if (!is.data.frame(fit_regions_df)) {
        fit_regions_df <-
          read.table(
            file = fit_regions_df,
            stringsAsFactors = FALSE,
            col.names = c("chr", "start", "end", "counts")
          )
      } else {
        colnames(fit_regions_df) <- c("chr", "start", "end", "counts")
      }
      
      fit_regions_params_df <-
        pbsR:::getDistributionParametersWithOptim(fit_regions_df, theta)
      
      corr_pbs_obj = pbsR:::getProbabilityBeingSignal(bin_df, fit_regions_params_df, theta)
      corr_pbs_obj$fit_regions_df = fit_regions_df
      corr_pbs_obj$fit_regions_lambda = fit_regions_params_df$lambda
      if (plot_pbs) {
        print(pbsR:::getPBSPlot(corr_pbs_obj))
      }
    }
    else{
      corr_pbs_obj = NULL
      if (plot_pbs) {
        print(pbsR:::getCorrectedPBSPlot(pbs_obj))
      }
    }
    return(list("pbs_obj" = pbs_obj, "corr_pbs_obj" = corr_pbs_obj))
  }
