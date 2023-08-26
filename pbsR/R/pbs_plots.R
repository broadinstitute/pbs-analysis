#' @importFrom ggplot2 ggplot aes geom_histogram geom_line ggtitle theme_bw scale_x_continuous after_stat
#' @importFrom rlang .data
getMultiplotHistograms <- function(max_dens, beta, k, lambda, working_df, bin_width) {
  sim_df <- data.frame(
    x = seq(0, max_dens, length.out = 500),
    y = dgamma(x = seq(0, max_dens, length.out = 500), shape = k, rate = beta)
  )
  # plot data
  plt <- ggplot(working_df, aes(x = .data$counts, y = after_stat(density))) +
    geom_histogram(position = "identity", binwidth = bin_width, alpha = 0.5) +
    geom_line(data = sim_df, aes(x = .data$x, y = .data$y * lambda)) +
    scale_x_continuous(limits = c(0, max_dens)) +
    ggtitle(label = paste("Bin counts histogram")) +
    theme_bw(base_size = 16)
  return(plt)
}

#' @importFrom ggplot2 ggplot aes geom_tile
#' @importFrom rlang .data
plotCVMHeatmap <- function(cvm_stat_df) {
  p = ggplot(cvm_stat_df) +
    geom_tile(aes(x = .data$beta, y = .data$k, fill = .data$omega2))
  
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_histogram geom_line theme_bw ggtitle scale_x_continuous labs scale_color_manual after_stat
#' @importFrom rlang .data
plotFittedDistribution <- function(working_df, empirical_density_df, fitted_density_df, max_dens, title_str, beta, k, lambda) {
  p <- ggplot(working_df, aes(x = .data$counts, y = after_stat(density))) +
    geom_histogram(binwidth = 0.05, fill = "grey66") +
    geom_line(data = empirical_density_df, aes(x = .data$x, y = .data$y, color = "lightgreen"), size = 1) +
    theme_bw(base_size = 16) +
    ggtitle(label = title_str, subtitle = paste0("beta, k, lambda: ", round(beta, 2), ", ", round(k, 2), ", ", round(lambda, 2))) +
    geom_line(data = fitted_density_df, aes(x = .data$x, y = lambda * .data$y, color = "salmon"), size = 1) +
    geom_line(data = empirical_density_df, aes(x = .data$x, y = .data$pbs, color = "lightblue"), size = 1) +
    scale_x_continuous(limits = c(0, max_dens)) +
    labs(x = "Counts", y = "Density") +
    scale_color_manual(
      labels = c("PBS", "Empirical", "Estimated"), name = "",
      values = c("orange", "blue", "mediumpurple")
    )
  return(p)
}

#find areas under curve between empirical and fitted curves, to the left and right of theta.
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_line geom_vline labs 
#' @importFrom rlang .data
plotAUC <- function(working_df, params_df, xlim = 20, theta = 0.5) { 
  
  densityPts <- density(working_df$counts, n = 1024, to = xlim)
  
  fittedPts <- dgamma(x = densityPts$x, shape = params_df$k, rate = params_df$beta)
  fittedPts <- as.data.frame(fittedPts)
  names(fittedPts) <- c("yFitted")
  
  lambda <- params_df$lambda
  
  difference_df <- data.frame(xValues = densityPts$x, yEmpirical = densityPts$y, yFitted = lambda*fittedPts)
  areas <- data.frame(LH = double(1), RH = double(1))
  
  #Right-Hand Area
  filtered_difference <- difference_df %>% dplyr::filter(xValues > quantile(working_df$counts, theta) ) %>% dplyr::filter(yEmpirical < yFitted)
  filtered_difference$difference <- (filtered_difference$yFitted - filtered_difference$yEmpirical)
  rhArea <- sum(filtered_difference$difference)
  
  #Left-Hand Area
  filtered_difference <- difference_df %>% dplyr::filter(xValues < quantile(working_df$counts, theta) )
  filtered_difference$difference <- abs(filtered_difference$yFitted - filtered_difference$yEmpirical)
  lhArea <- sum(filtered_difference$difference)
  
  diffPlot <- ggplot(data = difference_df, aes(x = xValues, y = yEmpirical, col = "Empirical")) + geom_line() + geom_line(data = difference_df, aes(y = yFitted, col = "Fitted"))
  p = diffPlot + geom_vline(xintercept = quantile(working_df$counts, theta), linetype = "dotted") + labs(x = "Counts", y = "Density")
  return(p)
}
