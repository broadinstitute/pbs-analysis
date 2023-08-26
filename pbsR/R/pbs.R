#' @export
PBS <- function(bin_df, theta = 0.5) {
  if (!is.data.frame(bin_df)) {
    bin_df <- read.table(file = bin_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(bin_df) <- c("chr", "start", "end", "counts")
  }
  params_df <- pbsR:::getDistributionParametersWithOptim(bin_df, theta = theta)
  return(pbsR:::getProbabilityBeingSignal(bin_df, params_df, return_plt = TRUE))
}

#' @export
iterPBS <- function(bin_df, a_comp_df, b_comp_df, theta = 0.5) {
  if (!is.data.frame(bin_df)) {
    bin_df <- read.table(file = bin_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(bin_df) <- c("chr", "start", "end", "counts")
  }
  pr1 <- pbsR:::getDistributionParametersWithOptim(bin_df, theta)

  p1 <- pbsR:::getMultiplotHistograms(20, pr1$beta, pr1$k, pr1$lambda, bin_df, 0.05) + ggtitle(label = "Bin Counts histogram - All", subtitle = paste0("beta, k, lambda: ", round(pr1[1, 1], 2), ", ", round(pr1[1, 2], 2), ", ", round(pr1[1, 3], 2)))

  if (!is.data.frame(a_comp_df)) {
    a_comp_df <- read.table(file = a_comp_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(a_comp_df) <- c("chr", "start", "end", "counts")
  }

  pr2 <- pbsR:::getDistributionParametersWithOptim(a_comp_df, lambda_range = c(0.5, 1), length_out = 10, max_dens = 20, theta)

  p2 <- pbsR:::getMultiplotHistograms(20, pr2$beta, pr2$k, pr2$lambda, a_comp_df, 0.05) + ggtitle(label = "Bin Counts histogram - A comp", subtitle = paste0("beta, k, lambda: ", round(pr2[1, 1], 2), ", ", round(pr2[1, 2], 2), ", ", round(pr2[1, 3], 2)))

  if (!is.data.frame(b_comp_df)) {
    b_comp_df <- read.table(file = b_comp_df, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "counts"))
  } else {
    colnames(b_comp_df) <- c("chr", "start", "end", "counts")
  }
  pr3 <- pbsR:::getDistributionParametersWithOptim(b_comp_df, lambda_range = c(0.5, 1), length_out = 10, max_dens = 20, theta)

  p3 <- pbsR:::getMultiplotHistograms(20, pr3$beta, pr3$k, pr3$lambda, b_comp_df, 0.05) + ggtitle(label = "Bin Counts histogram - B comp", subtitle = paste0("beta, k, lambda: ", round(pr3[1, 1], 2), ", ", round(pr3[1, 2], 2), ", ", round(pr3[1, 3], 2)))

  fit_per_compartment <- patchwork::wrap_plots(p1, p2, p3, ncol = 1)

  return(list(
    "pbs_obj" = pbsR:::getProbabilityBeingSignal(bin_df, pr1, return_plt = TRUE),
    "corr_pbs_obj" = pbsR:::getProbabilityBeingSignal(bin_df, pr2, return_plt = TRUE),
    "fit_per_compartment" = fit_per_compartment
  ))
}
