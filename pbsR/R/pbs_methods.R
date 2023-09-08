#' @export
getProbabilityBeingSignal <- function(bin_df,
                                      params_df,
                                      theta,
                                      max_dens = 20) {
  working_df <- bin_df[bin_df$counts > 0,]
  threshold <-
    quantile(x = working_df$counts[working_df$chr %in% paste0("chr", 1:22)], probs = 0.9999)
  APP_density <-
    density(x = working_df$counts[working_df$counts < threshold],
            adjust = 1,
            n = 10000)
  empirical_density_df <-
    data.frame(x = APP_density$x, y = APP_density$y)
  
  # get beta, k, lambda from params_df
  beta = params_df$beta
  k = params_df$k
  
  p_values <-
    pgamma(
      q = working_df$counts,
      rate = beta,
      shape = k,
      lower.tail = FALSE
    )
  ratio <- mean(p_values > (1 - theta))
  lambda <-  ratio / theta
  
  #lambda = params_df$lambda
  
  print(paste0("########## Assigning PBS score to each bin #################"))
  
  # calculate density for background
  fitted_density_df <-
    data.frame(x = APP_density$x,
               y = dgamma(
                 x = APP_density$x,
                 shape = k,
                 rate = beta
               ))
  empirical_density_df$pbs <-
    (empirical_density_df$y - lambda * fitted_density_df$y) / empirical_density_df$y
  
  # set values of pbs to be zero where it looks unstable
  first_zero_idx <-
    max(which(
      empirical_density_df$pbs <= 0 &
        is.finite(empirical_density_df$pbs)
    ))
  empirical_density_df$pbs[1:first_zero_idx] <- 0
  
  # set values of pbs to equal one where they're greater than first_zero_idx and are not finite
  empirical_density_df$pbs[empirical_density_df$pbs > first_zero_idx &
                             !is.finite(empirical_density_df$pbs)] <- 1
  
  # add some points to empirical_density_df to represent highest value of counts
  empirical_density_df <-
    rbind(empirical_density_df, data.frame(
      x = max(working_df$counts, na.rm = TRUE),
      y = 0,
      pbs = 1
    ))
  
  # need to approximate values in working_df based on empirical_density_pbs
  estimated_pbs <-
    stats::approx(x = empirical_density_df$x,
                  y = empirical_density_df$pbs,
                  xout = bin_df$counts)
  
  # some NAs in estimated_pbs: set these to be zero
  bin_df$pbs <- estimated_pbs[[2]]
  bin_df$pbs[is.na(bin_df$pbs)] <- 0
  
  #if (is.null(title_str)) {
  # title_str <- paste("PBS")
  #}
  return(
    list(
      "bin_df" = bin_df,
      "empirical_density_df" = empirical_density_df,
      "fitted_density_df" = fitted_density_df,
      "beta" = beta,
      "k" = k,
      "lambda" = lambda,
      "RH_area" = params_df$RH_area
    )
  )
}
