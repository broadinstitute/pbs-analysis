#' @export
getProbabilityBeingSignal <- function(bin_df, params_df,
                                      plot_data = FALSE, title_str = NULL,
                                      return_plt = FALSE, plot_distribution_fit_only = FALSE, max_dens = 10) {
  working_df <- bin_df[bin_df$counts > 0, ]
  threshold <- quantile(x = working_df$counts[working_df$chr %in% paste0("chr", 1:22)], probs = 0.9999)
  APP_density <- density(x = working_df$counts[working_df$counts < threshold], adjust = 1, n = 10000)
  empirical_density_df <- data.frame(x = APP_density$x, y = APP_density$y)
  # get beta and k from params_df
  beta <- params_df$beta
  k <- params_df$k
  # calculate lambda based on ratio between left and right sides of curve (i.e. what fraction of points lie below 0.5 vs above)
  p_values <- pgamma(q = working_df$counts, shape = params_df$k, rate = params_df$beta, lower.tail = FALSE)
  ratio <- mean(p_values > 0.5)
  lambda <- 2 * ratio
  # calculate density for background
  fitted_density_df <- data.frame(x = APP_density$x, y = dgamma(x = APP_density$x, shape = k, rate = beta))
  empirical_density_df$pbs <- (empirical_density_df$y - lambda * fitted_density_df$y) / empirical_density_df$y
  # set values of pbs to be zero where it looks unstable
  first_zero_idx <- max(which(empirical_density_df$pbs <= 0 & is.finite(empirical_density_df$pbs)))
  empirical_density_df$pbs[1:first_zero_idx] <- 0
  # set values of pbs to equal one where they're greater than first_zero_idx and are not finite
  empirical_density_df$pbs[empirical_density_df$pbs > first_zero_idx & !is.finite(empirical_density_df$pbs)] <- 1
  # add some points to empirical_density_df to represent highest value of counts
  empirical_density_df <- rbind(empirical_density_df, data.frame(x = max(working_df$counts, na.rm = TRUE), y = 0, pbs = 1))
  # need to approximate values in working_df based on empirical_density_pbs
  estimated_pbs <- stats::approx(x = empirical_density_df$x, y = empirical_density_df$pbs, xout = bin_df$counts)
  # some NAs in estimated_pbs: set these to be zero
  bin_df$pbs <- estimated_pbs[[2]]
  bin_df$pbs[is.na(bin_df$pbs)] <- 0

  if (is.null(title_str)) {
    title_str <- paste("PBS")
  }
  if (plot_distribution_fit_only) {
    median_countss <- median(working_df$counts)
  } else {
    plt1 <- pbsR:::plotFittedDistribution(working_df, empirical_density_df, fitted_density_df, max_dens, title_str, beta, k, lambda)
  }
  if (plot_data) {
    print(plt1)
  }
  if (return_plt) {
    return(list("bin_df" = bin_df, "pbs_plt" = plt1))
  } else {
    return(bin_df)
  }
}
