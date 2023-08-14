#' @export
PBS <- function(bin_df) {
  if (!is.data.frame(bin_df)) {
    bin_df <- read.table(file = bin_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(bin_df) <- c("chr", "start", "end", "counts")
  }
  params_df <- getDistributionParametersWithOptim(bin_df)
  return(getProbabilityBeingSignal(bin_df, params_df, return_plt = TRUE))
}

#' @export
writePBS <- function(pbs_obj, output_filename) {
  write.table(file = output_filename, x = pbs_obj$bin_df[, c("chr", "start", "end", "pbs")], quote = F, row.names = F, sep = "\t")
}

#' @export
iterPBS <- function(bin_df, a_comp_df, b_comp_df) {
  if (!is.data.frame(bin_df)) {
    bin_df <- read.table(input = bin_df, data.table = FALSE, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(bin_df) <- c("chr", "start", "end", "counts")
  }
  pr1 <- getDistributionParametersWithOptim(bin_df)
  p1 <- getMultiplotHistograms(20, pr1$beta, pr1$k, pr1$lambda, df, 0.05) + ggtitle(label = "Bin Counts histogram - All", subtitle = paste0("beta, k, lambda: ", round(pr1[1, 1], 2), ", ", round(pr1[1, 2], 2), ", ", round(pr1[1, 3], 2)))

  if (!is.data.frame(a_comp_df)) {
    a_comp_df <- read.table(input = a_comp_df, data.table = FALSE, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(a_comp_df) <- c("chr", "start", "end", "counts")
  }

  pr2 <- getDistributionParametersWithOptim(a_comp_df)
  p2 <- getMultiplotHistograms(20, pr2$beta, pr2$k, pr2$lambda, adf, 0.05) + ggtitle(label = "Bin Counts histogram - A comp", subtitle = paste0("beta, k, lambda: ", round(pr2[1, 1], 2), ", ", round(pr2[1, 2], 2), ", ", round(pr2[1, 3], 2)))

  if (!is.data.frame(b_comp_df)) {
    b_comp_df <- read.table(input = b_comp_df, data.table = FALSE, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(b_comp_df) <- c("chr", "start", "end", "counts")
  }
  pr3 <- getDistributionParametersWithOptim(b_comp_df)
  p3 <- getMultiplotHistograms(20, pr3$beta, pr3$k, pr3$lambda, bdf, 0.05) + ggtitle(label = "Bin Counts histogram - B comp", subtitle = paste0("beta, k, lambda: ", round(pr3[1, 1], 2), ", ", round(pr3[1, 2], 2), ", ", round(pr3[1, 3], 2)))

  fit_per_compartment <- p1 + p2 + p3 + plot_layout(ncol = 1)

  pbs_obj <- getProbabilityBeingSignal(bin_df, pr1, return_plt = TRUE)
  corr_pbs_obj <- getProbabilityBeingSignal(bin_df, pr2, return_plt = TRUE)

  return(list(
    "bin_df" = pbs_obj$bin_df,
    "pbs_plt" = pbs_obj$plt1,
    "corr_bin_df" = corr_pbs_obj$bin_df,
    "corr_pbs_plt" = corr_pbs_obj$plt1
  ))
}

#' @export
getPBS <- function(pbs_obj) {
  return(pbs_obj$bin_df)
}

#' @export
plotPBS <- function(pbs_obj) {
  print(pbs_obj$pbs_plt)
}

#' @export
getCorrectedPBS <- function(pbs_obj) {
  return(pbs_obj$corr_bin_df)
}

#' @export
plotCorrectedPBS <- function(pbs_obj) {
  print(pbs_obj$corr_pbs_plt)
}
