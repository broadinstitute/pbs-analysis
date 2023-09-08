#' @export
getPBS <- function(pbs_obj) {
  return(pbs_obj$pbs_obj$bin_df)
}

#' @export
getCorrectedPBS <- function(pbs_obj) {
  return(pbs_obj$corr_pbs_obj$bin_df)
}

#' @export
getPBSPlot <- function(pbs_obj,
                       max_dens = 20,
                       title_str = "PBS") {
  plot_obj = pbs_obj
  if(!is.null(pbs_obj$pbs_obj)){
    plot_obj = pbs_obj$pbs_obj
  }
  pbsR:::plotFittedDistribution(
    plot_obj$bin_df,
    plot_obj$empirical_density_df,
    plot_obj$fitted_density_df,
    max_dens = max_dens,
    title_str = title_str,
    plot_obj$beta,
    plot_obj$k,
    plot_obj$lambda
  )
}

#' @export
getCorrectedPBSPlot <-
  function(pbs_obj,
           max_dens = 20,
           title_str = "PBS") {
    plot_obj = pbs_obj
    if(!is.null(pbs_obj$corr_pbs_obj)){
      plot_obj = pbs_obj$corr_pbs_obj
    }
    pbsR:::plotFittedDistribution(
      plot_obj$bin_df,
      plot_obj$empirical_density_df,
      plot_obj$fitted_density_df,
      max_dens = max_dens,
      title_str = title_str,
      plot_obj$beta,
      plot_obj$k,
      plot_obj$lambda
    )
  }

#' @export
getCompartmentFits <- function(pbs_obj) {
  pbsR:::plotMultiplotHistograms(pbs_obj)
}

#' @export
writePBS <- function(pbs_obj, output_filename) {
  write.table(
    file = output_filename,
    x = pbs_obj$pbs_obj$bin_df,
    quote = F,
    row.names = F,
    sep = "\t"
  )
}

#' @export
writeCorrectedPBS <- function(pbs_obj, output_filename) {
  write.table(
    file = output_filename,
    x = pbs_obj$corr_pbs_obj$bin_df,
    quote = F,
    row.names = F,
    sep = "\t"
  )
}
