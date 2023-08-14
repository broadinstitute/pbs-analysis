#' @export
PBS <- function(bin_df) {
  if (!is.data.frame(bin_df)) {
    bin_df <- read.table(file = bin_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(bin_df) <- c("chr", "start", "end", "counts")
  }
  params_df <- pbsR:::getDistributionParametersWithOptim(bin_df)
  return(pbsR:::getProbabilityBeingSignal(bin_df, params_df, return_plt = TRUE))
  #pbs_obj = pbsR:::getProbabilityBeingSignal(bin_df, params_df, return_plt = TRUE)
  #return(list("bin_df" = pbs_obj$bin_df, "pbs_plt" = pbs_obj$plt1))
}

#' @export
writePBS <- function(pbs_obj, output_filename) {
  write.table(file = output_filename, x = pbs_obj$bin_df[, c("chr", "start", "end", "pbs")], quote = F, row.names = F, sep = "\t")
}

#' @export
iterPBS <- function(bin_df, a_comp_df, b_comp_df) {
  if (!is.data.frame(bin_df)) {
    bin_df <- read.table(file = bin_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(bin_df) <- c("chr", "start", "end", "counts")
  }
  pr1 <- pbsR:::getDistributionParametersWithOptim(bin_df)
  p1 <- pbsR:::getMultiplotHistograms(20, pr1$beta, pr1$k, pr1$lambda, bin_df, 0.05) + ggtitle(label = "Bin Counts histogram - All", subtitle = paste0("beta, k, lambda: ", round(pr1[1, 1], 2), ", ", round(pr1[1, 2], 2), ", ", round(pr1[1, 3], 2)))

  if (!is.data.frame(a_comp_df)) {
    a_comp_df <- read.table(file = a_comp_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(a_comp_df) <- c("chr", "start", "end", "counts")
  }

  pr2 <- pbsR:::getDistributionParametersWithOptim(a_comp_df)
  p2 <- pbsR:::getMultiplotHistograms(20, pr2$beta, pr2$k, pr2$lambda, a_comp_df, 0.05) + ggtitle(label = "Bin Counts histogram - A comp", subtitle = paste0("beta, k, lambda: ", round(pr2[1, 1], 2), ", ", round(pr2[1, 2], 2), ", ", round(pr2[1, 3], 2)))

  if (!is.data.frame(b_comp_df)) {
    b_comp_df <- read.table(file = b_comp_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(b_comp_df) <- c("chr", "start", "end", "counts")
  }
  pr3 <- pbsR:::getDistributionParametersWithOptim(b_comp_df)
  p3 <- pbsR:::getMultiplotHistograms(20, pr3$beta, pr3$k, pr3$lambda, b_comp_df, 0.05) + ggtitle(label = "Bin Counts histogram - B comp", subtitle = paste0("beta, k, lambda: ", round(pr3[1, 1], 2), ", ", round(pr3[1, 2], 2), ", ", round(pr3[1, 3], 2)))

  #fit_per_compartment <- gridExtra::grid.arrange(p1,p2,p3,ncol=1)
  fit_per_compartment = patchwork::wrap_plots(p1,p2,p3, ncol = 1)

  #pbs_obj <- pbsR:::getProbabilityBeingSignal(bin_df, pr1, return_plt = TRUE)
  #corr_pbs_obj <- pbsR:::getProbabilityBeingSignal(bin_df, pr2, return_plt = TRUE)

  return(list(
    "pbs_obj" = pbsR:::getProbabilityBeingSignal(bin_df, pr1, return_plt = TRUE),
    "corr_pbs_obj" = pbsR:::getProbabilityBeingSignal(bin_df, pr2, return_plt = TRUE),
    "fit_per_compartment" = fit_per_compartment
  ))
}

#' @export
getPBS <- function(pbs_obj) {
  if("bin_df" %in% names(pbs_obj)){
    return(pbs_obj$bin_df)
  }
  else{
    getPBS(pbs_obj$pbs_obj)
  }
}

#' @export
plotPBS <- function(pbs_obj) {
  if("pbs_plt" %in% names(pbs_obj)){
    print(pbs_obj$pbs_plt)
  }
  else{
    plotPBS(pbs_obj$pbs_obj)
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
