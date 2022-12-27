#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
GetProbabilityBeingSignal <- function(bin_df_filename, params_df_filename, pbs_filename = NULL,
                                 plot_data = FALSE, title_str = NULL, return_bin_df = TRUE,
                                 return_plt = FALSE, plot_distribution_fit_only = FALSE, max_dens = 10){
  if(!is.data.frame(bin_df_filename)){
    bin_df <- fread(input = bin_df_filename, data.table = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'start', 'end', 'counts'))
  } else{
    bin_df <- bin_df_filename
  }
  if(!is.data.frame(params_df_filename)){
    params_df <- fread(input = params_df_filename, data.table = FALSE, stringsAsFactors = FALSE)
  } else{
    params_df <- params_df_filename
  }
  working_df <- bin_df[bin_df$counts > 0,]
  threshold <- quantile(x = working_df$counts[working_df$chr %in% paste0('chr', 1:22)], probs = 0.9999)
  APP_density <- density(x = working_df$counts[working_df$counts < threshold], adjust = 1, n = 10000)
  empirical_density_df <- data.frame(x = APP_density$x, y = APP_density$y)
  # get beta and k from params_df
  beta <- params_df$beta
  k <- params_df$k
  # calculate lambda based on ratio between left and right sides of curve (i.e. what fraction of points lie below 0.5 vs above)
  p_values <- pgamma(q = working_df$counts, shape = params_df$k, rate = params_df$beta, lower.tail = FALSE)
  ratio <- mean(p_values > 0.5)
  lambda <- 2*ratio
  # calculate density for background
  fitted_density_df <- data.frame(x = APP_density$x, y = dgamma(x = APP_density$x, shape = k, rate = beta))
  empirical_density_df$pbs <- (empirical_density_df$y - lambda*fitted_density_df$y)/empirical_density_df$y
  # set values of pbs to be zero where it looks unstable
  first_zero_idx <- max(which(empirical_density_df$pbs <= 0 & is.finite(empirical_density_df$pbs)))
  empirical_density_df$pbs[1:first_zero_idx] <- 0
  # set values of pbs to equal one where they're greater than first_zero_idx and are not finite
  empirical_density_df$pbs[empirical_density_df$pbs > first_zero_idx & !is.finite(empirical_density_df$pbs)] <- 1
  # add some points to empirical_density_df to represent highest value of counts
  empirical_density_df <- rbind(empirical_density_df, data.frame(x = max(working_df$counts, na.rm = TRUE), y = 0, pbs = 1))
  # need to approximate values in working_df based on empirical_density_pbs
  estimated_pbs <- approx(x = empirical_density_df$x, y = empirical_density_df$pbs, xout = bin_df$counts)
  # some NAs in estimated_pbs: set these to be zero
  bin_df$pbs <- estimated_pbs[[2]]
  bin_df$pbs[is.na(bin_df$pbs)] <- 0
  if(!is.null(pbs_filename)){
    write.table(x = bin_df[,c('chr', 'start', 'end', 'pbs')], file = pbs_filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  if(is.null(title_str)){
    title_str <- paste('PBS', basename(bin_df_filename))
  }
  if(plot_distribution_fit_only){
    median_countss <- median(working_df$counts)
    plt1 <- ggplot(working_df, aes(x = counts, y = ..density..)) + geom_histogram(binwidth = 0.01, fill = 'grey66') +
      geom_line(data = empirical_density_df, aes(x = x, y = y, color = 'lightgreen'), size = 1) +
      theme_bw(base_size = 24) + ggtitle(label = title_str) +
      geom_line(data = fitted_density_df, aes(x = x, y = lambda*y, color = 'salmon'), size = 1) +
      geom_vline(linetype = 'dotted', color = 'black', xintercept = median_countss) +
      scale_x_continuous(limits = c(0, max_dens)) + labs(x = 'Count', y = 'Density') +
      scale_color_manual(labels = c('Empirical', 'Estimated'), name = '',
                         values = c('blue', 'mediumpurple'))
  } else{
    plt1 <- ggplot(working_df, aes(x = counts, y = ..density..)) + geom_histogram(binwidth = 0.05, fill = 'grey66') +
      geom_line(data = empirical_density_df, aes(x = x, y = y, color = 'lightgreen'), size = 1) +
      theme_bw(base_size = 16) + ggtitle(label = title_str) +
      geom_line(data = fitted_density_df, aes(x = x, y = lambda*y, color = 'salmon'), size = 1) +
      geom_line(data = empirical_density_df, aes(x = x, y = pbs, color = 'lightblue'), size = 1) +
      scale_x_continuous(limits = c(0, max_dens)) + labs(x = 'Counts', y = 'Density') +
      scale_color_manual(labels = c('PBS', 'Empirical', 'Estimated'), name = '',
                         values = c('orange', 'blue', 'mediumpurple'))
  }
  if(plot_data){
    print(plt1)
  }
  if(return_bin_df & return_plt){
    return(list('bin_df' = bin_df, 'pbs_plt' = plt1))
  } else if(return_plt){
    return(plt1)
  } else if(return_bin_df){
    return(bin_df)
  } else{
    return(TRUE)
  }
}