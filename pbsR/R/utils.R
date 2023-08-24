#' @export
getPBS <- function(pbs_obj) {
  if ("bin_df" %in% names(pbs_obj)) {
    return(pbs_obj$bin_df)
  } else {
    getPBS(pbs_obj$pbs_obj)
  }
}

#' @export
plotPBS <- function(pbs_obj) {
  if ("pbs_plot" %in% names(pbs_obj)) {
    print(pbs_obj$pbs_plot)
  } else {
    plotPBS(pbs_obj$pbs_obj)
  }
}

#' @export
writePBS <- function(pbs_obj, output_filename) {
  if ("bin_df" %in% names(pbs_obj)) {
    write.table(file = output_filename, x = pbs_obj$bin_df, quote = F, row.names = F, sep = "\t")
  } else {
    writePBS(pbs_obj$corr_pbs_obj)
  }
}

#' @export
getCorrectedPBS <- function(pbs_obj) {
  getPBS(pbs_obj$corr_pbs_obj)
}

#' @export
plotCorrectedPBS <- function(pbs_obj) {
  plotPBS(pbs_obj$corr_pbs_obj)
}

#' @export
plotCompartmentFits <- function(pbs_obj) {
  print(pbs_obj$fit_per_compartment)
}
