#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(fBasics))
suppressPackageStartupMessages(library(goftest))

GetParameterCombinations <- function(beta_range = c(0.5, 1), k_range = c(18, 25), lambda_range = c(0.6, 0.9)){
  params_all <- expand.grid(beta = beta_range, k = k_range, lambda = lambda_range)
  # refine the grid a little
  params_all$params_num <- 1:nrow(params_all)
  return(params_all)
}

GetMultiplotHistograms <- function(max_dens, beta, k, lambda, working_df, bin_width){
  sim_df <- data.frame(x = seq(0, max_dens, length.out = 500),
                       y = dgamma(x = seq(0, max_dens, length.out = 500), shape = k, rate = beta))
  # plot data
  plt <- ggplot(working_df, aes(x = counts, y = ..density..)) +
    geom_histogram(position = 'identity', binwidth = bin_width, alpha = 0.5) +
    geom_line(data = sim_df, aes(x = x, y = y*lambda)) + scale_x_continuous(limits = c(0, max_dens)) +
    ggtitle(label = paste('Bin counts histogram')) + theme_bw(base_size = 16)
  print(plt)
}

cvm.test2 <- function(x, null="punif", ..., nullname) {
  xname <- deparse(substitute(x))
  nulltext <- deparse(substitute(null))
  if(is.character(null)) nulltext <- null
  if(missing(nullname) || is.null(nullname)) {
    reco <- recogniseCdf(nulltext)
    nullname <- if(!is.null(reco)) reco else
      paste("distribution", sQuote(nulltext))
  }
  stopifnot(is.numeric(x))
  x <- as.vector(x)
  n <- length(x)
  F0 <- if(is.function(null)) null else
    if(is.character(null)) get(null, mode="function") else
      stop("Argument 'null' should be a function, or the name of a function")
  U <- F0(x, ...)
  if(any(is.nan(U)) | any(is.na(U))){
    browser()
  }
  if(any(U < 0 | U > 1))
    U[U < 0] = 0; U[U > 1] = 1
  # print('Warning: coercing U to between 0 and 1')
  #     stop("null distribution function returned values outside [0,1]")
  U <- sort(U)
  k <- seq_len(n)
  omega2 <- 1/(12 * n) + sum((1-Heaviside(k,0.5*n))*(U - (2 * k - 1)/(2 * n))^2)
  PVAL <- pCvM(omega2, n=n, lower.tail=FALSE)
  names(omega2) <- "omega2"
  METHOD <- c("Cramer-von Mises test of goodness-of-fit",
              paste("Null hypothesis:", nullname))
  extras <- list(...)
  parnames <- intersect(names(extras), names(formals(F0)))
  if(length(parnames) > 0) {
    pars <- extras[parnames]
    pard <- character(0)
    for(i in seq_along(parnames))
      pard[i] <- paste(parnames[i], "=", paste(pars[[i]], collapse=" "))
    pard <- paste("with",
                  ngettext(length(pard), "parameter", "parameters"),
                  "  ",
                  paste(pard, collapse=", "))
    METHOD <- c(METHOD, pard)
  }
  out <- list(statistic = omega2,
              p.value = PVAL,
              method = METHOD,
              data.name = xname)
  class(out) <- "htest"
  return(out)
}

pgamma_null <- function(q, params_df, weight_value){
  pgamma_mix <- params_df$lambda*pgamma(q = q, shape = params_df$k, rate = params_df$beta) +
    (1 - params_df$lambda)*punif(q = q, min = weight_value, max = 1.01*weight_value)
  return(pgamma_mix)
}

# try calculating cost using Cramer Von Mises
getCVMDistance <- function(working_df, params_df, weight_value){

  if(!"counts" %in% colnames(working_df)) {
    names(working_df) <- c("Chr", "StartIdx", "counts")
  }
  cvm_stat <- cvm.test2(x = working_df$counts, null = 'pgamma_null', 'params_df' = params_df,
                        'weight_value' = weight_value)
  return(cvm_stat$statistic)
}

GetDistributionParameters <- function(working_df, lambda_range = seq(from = 0.6, to = 1, by = 0.05),
                                      fix_weight = TRUE, weight_value = 500, use_log = FALSE, length_out = 100,
                                      max_range_multiple = 300, plot_data = FALSE){
  working_df <- working_df %>% dplyr::filter(counts > 0)
  if(use_log){
    working_df$log_count <- log(working_df$counts)
    working_df$counts <- working_df$log_count
    bin_width <- 0.1; max_dens <- 10
  } else{
    bin_width <- 10; max_dens <- 500
  }
  mean_init <- mean(working_df$counts); var_init <- var(working_df$counts)
  beta_init <- mean_init/var_init; k_init <- mean_init*beta_init
  params_all <- GetParameterCombinations(beta_range = seq(beta_init, max_range_multiple*beta_init, length.out = length_out),
                                         k_range = seq(k_init, max_range_multiple*k_init, length.out = length_out),
                                         lambda_range = lambda_range)
  # set max and min vals for weight at far limit of distribution
  if(!fix_weight){
    weight_value <- 0.99*max(working_df$counts)
  }
  cvm_stat_list <- ddply(.data = params_all, .variables = 'params_num', .fun = getCVMDistance, 'working_df' = working_df,
                         'weight_value' = weight_value)
  # make a heat map of cvm_stat_list
  cvm_stat_list <- inner_join(x = cvm_stat_list, y = params_all, by = "params_num")
  if(plot_data){
    ggplot(cvm_stat_list, aes(x = beta, y = k)) + geom_tile(aes(fill = omega2))
  }
  # what's the parameter value at the minimum cost?
  optim_params <- params_all[which.min(cvm_stat_list$omega2),]
  optim_params$Name <- working_df$Name[1]
  return(optim_params)
}

# get initial values from GetDistributionParameters with coarse-grained grid
GetDistributionParametersWithOptim <- function(working_df, lambda_range = seq(0.6, 0.9, by = 0.1), length_out = 3,
                                               fix_weight = TRUE, weight_value = 500, plot_data = FALSE, plot_terra = FALSE,
                                               bin_width = 0.05, max_dens = 10){
  counts_col <- grepl(pattern = 'count', x = names(working_df), ignore.case = TRUE)
  if(sum(counts_col) != 1){
    stop('Working_df must have exactly one column with a name similar to "counts".')
  }
  names(working_df)[counts_col] <- 'counts'
  working_df <- working_df %>% dplyr::filter(counts > 0)
  params_init <- GetDistributionParameters(working_df = working_df, lambda_range = lambda_range, length_out = length_out)
  # set max and min vals for weight at far limit of distribution
  beta_init <- params_init$beta; k_init <- params_init$k; lambda_init <- params_init$lambda
  if(!fix_weight){
    weight_value <- 0.99*max(working_df$counts)
  }
  # use optim to calculate parameters
  optim_params <- optim(par = c(beta_init, k_init, lambda_init),
                        fn = function(par){params_df <- data.frame(beta = par[1], k = par[2], lambda = par[3])
                        getCVMDistance(params_df = params_df,
                                       weight_value = weight_value, working_df = working_df)},
                        method = 'L-BFGS-B', lower = c(1e-5, 1e-5, 1e-5), upper = c(Inf, Inf, 1))
  # calculate lambda based on  p values
  lambda <- 2*mean(pgamma(q = working_df$counts, shape = optim_params$par[2],
                          rate = optim_params$par[1], lower.tail = FALSE) > 0.5)
  # plot data
  if(plot_data){
    GetMultiplotHistograms(max_dens = max_dens, beta = optim_params$par[1], k = optim_params$par[2],
                           lambda = lambda, working_df = working_df, bin_width = bin_width)
  }
  # plot QC figure
  if(plot_terra){
    makeQCPlot(working_df=working_df, beta = optim_params$par[1], k = optim_params$par[2])
  }
  # get fit quality metric
  params_df <- data.frame(beta = optim_params$par[1], k = optim_params$par[2], lambda = lambda)
  fit_quality <- findAreasWithPlot(working_df = working_df, params_df = params_df, xlim = max_dens, printGraph = FALSE)
  params_df$resid_area <- fit_quality$RH
  return(params_df)
}

# get fit quality metric (from Sreshtaa Rajesh)
findAreasWithPlot <- function(working_df, params_df, xlim = 10, filetitle = NULL, printGraph = FALSE) { 
  densityPts <- density(working_df$counts, n = 1024, to = xlim)
  fittedPts <- dgamma(x = densityPts$x, shape = params_df$k, rate = params_df$beta)
  fittedPts <- as.data.frame(fittedPts)
  names(fittedPts) <- c("yFitted")
  lambda <- params_df$lambda
  difference_df <- data.frame(xValues = densityPts$x, yEmpirical = densityPts$y, yFitted = lambda*fittedPts)
  areas <- data.frame(RH = double(1))
  
  #Right-Hand Area
  filtered_difference <- difference_df %>% dplyr::filter(xValues > median(working_df$counts)) %>% 
                            dplyr::filter(yEmpirical < yFitted)
  filtered_difference$difference <- (filtered_difference$yFitted - filtered_difference$yEmpirical)
  rhArea <- sum(filtered_difference$difference)*(densityPts$x[2] - densityPts$x[1])
  areas[1,] <- rhArea
  return(areas)
}

# Produce QC plot (for ENCODE Workshop)
makeQCPlot <- function(working_df, beta, k){
  working_df <- working_df[working_df$counts > 0, , drop=F]
  threshold <- quantile(x = working_df$counts, probs = 0.9999)
  APP_density <- density(x = working_df$counts[working_df$counts < threshold], adjust = 1, n = 10000)
  empirical_density_df <- data.frame(x = APP_density$x, y = APP_density$y)

  # calculate lambda based on ratio between left and right sides of curve (i.e. what fraction of points lie below 0.5 vs above)
  p_values <- pgamma(q = working_df$counts, shape = k, rate = beta, lower.tail = FALSE)
  ratio <- mean(p_values > 0.5)
  lambda <- 2*ratio
  # calculate density for background
  fitted_density_df <- data.frame(x = APP_density$x, y = dgamma(x = APP_density$x, shape = k, rate = beta))
  # plot fitted density and empirical density
  empirical_density_df$metric <- (empirical_density_df$y - lambda*fitted_density_df$y)/empirical_density_df$y
  # set values of metric to be zero where it looks unstable
  first_zero_idx <- max(which(empirical_density_df$metric <= 0 & is.finite(empirical_density_df$metric)))
  empirical_density_df$metric[1:first_zero_idx] <- 0
  # set values of metric to equal one where they're greater than first_zero_idx and are not finite
  empirical_density_df$metric[empirical_density_df$metric > first_zero_idx & !is.finite(empirical_density_df$metric)] <- 1

  # Plot
  plt <- ggplot(working_df, aes(x = counts, y = ..density..)) + geom_histogram(binwidth = 0.05) + 
    geom_line(data = empirical_density_df, aes(x = x, y = y, color = 'darkred'), size = 1) + 
    theme_bw() + theme(text = element_text(size = 12)) + ggtitle(label = 'Distribution Fit') +
    geom_line(data = fitted_density_df, aes(x = x, y = lambda*y, color = 'darkblue'), size = 1) +
    geom_line(data = empirical_density_df, aes(x = x, y = metric, color = 'green'), size = 1) +
    scale_x_continuous(limits = c(0, 10)) + labs(x = 'Count', y = 'Density') +
    scale_color_manual(labels = c('Estimated', 'Empirical', 'Signal \nprobability'), name = '', 
                       values = c('salmon', 'turquoise', 'mediumpurple'))
  pdf('fit.pdf')
  print(plt)
  dev.off()
}