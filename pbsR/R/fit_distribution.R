getParameterCombinations <- function(beta_range = c(0.5, 1), k_range = c(18, 25), lambda_range = c(0.6, 0.9)) {
  params_all <- expand.grid(beta = beta_range, k = k_range, lambda = lambda_range)
  params_all$params_num <- 1:nrow(params_all)
  return(params_all)
}

cvm.test2 <- function(x, null = "punif", ..., nullname) {
  xname <- deparse(substitute(x))
  nulltext <- deparse(substitute(null))
  if (is.character(null)) nulltext <- null
  if (missing(nullname) || is.null(nullname)) {
    reco <- goftest::recogniseCdf(nulltext)
    nullname <- if (!is.null(reco)) {
      reco
    } else {
      paste("distribution", sQuote(nulltext))
    }
  }
  stopifnot(is.numeric(x))
  x <- as.vector(x)
  n <- length(x)
  F0 <- if (is.function(null)) {
    null
  } else if (is.character(null)) {
    get(null, mode = "function")
  } else {
    stop("Argument 'null' should be a function, or the name of a function")
  }
  U <- F0(x, ...)
  if (any(is.nan(U)) | any(is.na(U))) {
    browser()
  }
  if (any(U < 0 | U > 1)) {
    U[U < 0] <- 0
  }
  U[U > 1] <- 1
  # print('Warning: coercing U to between 0 and 1')
  #     stop("null distribution function returned values outside [0,1]")
  U <- sort(U)
  k <- seq_len(n)
  omega2 <- 1 / (12 * n) + sum((1 - fBasics::Heaviside(k, 0.5 * n)) * (U - (2 * k - 1) / (2 * n))^2)
  PVAL <- goftest::pCvM(omega2, n = n, lower.tail = FALSE)
  names(omega2) <- "omega2"
  METHOD <- c(
    "Cramer-von Mises test of goodness-of-fit",
    paste("Null hypothesis:", nullname)
  )
  extras <- list(...)
  parnames <- intersect(names(extras), names(formals(F0)))
  if (length(parnames) > 0) {
    pars <- extras[parnames]
    pard <- character(0)
    for (i in seq_along(parnames)) {
      pard[i] <- paste(parnames[i], "=", paste(pars[[i]], collapse = " "))
    }
    pard <- paste(
      "with",
      ngettext(length(pard), "parameter", "parameters"),
      "  ",
      paste(pard, collapse = ", ")
    )
    METHOD <- c(METHOD, pard)
  }
  out <- list(
    statistic = omega2,
    p.value = PVAL,
    method = METHOD,
    data.name = xname
  )
  class(out) <- "htest"
  return(out)
}

pgamma.null <- function(q, params_df, weight_value) {
  pgamma_mix <- params_df$lambda * stats::pgamma(q = q, shape = params_df$k, rate = params_df$beta) +
    (1 - params_df$lambda) * stats::punif(q = q, min = weight_value, max = 1.01 * weight_value)
  return(pgamma_mix)
}

# try calculating cost using Cramer Von Mises
getCVMDistance <- function(working_df, params_df, weight_value) {
  cvm_stat <- cvm.test2(x = working_df$counts, null = "pgamma.null", "params_df" = params_df, "weight_value" = weight_value)
  return(cvm_stat$statistic)
}

# run Cramer Von Mises
#' @importFrom magrittr %>%
runCVM <- function(working_df, params_df) {
  working_df = working_df %>% dplyr::filter(counts > 0)
  weight_value = 0.99*max(working_df$counts)
  cvm_stat <- cvm.test2(x = working_df$counts, null = "pgamma.null", "params_df" = params_df, "weight_value" = weight_value)
  print(c(cvm_stat$statistic, cvm_stat$p.value))
}

#' @export
#' @importFrom magrittr %>%
getDistributionParameters <- function(working_df, lambda_range = seq(from = 0.4, to = 0.8, by = 0.05),
                                      fix_weight = TRUE, weight_value = 500, use_log = FALSE, length_out = 100,
                                      max_range_multiple = 300, plot_data = FALSE) {
  
  working_df <- working_df %>% dplyr::filter(counts > 0)
  
  if (use_log) {
    working_df$log_count <- log(working_df$counts)
    working_df$counts <- working_df$log_count
    bin_width <- 0.1
    max_dens <- 10
  } else {
    bin_width <- 10
    max_dens <- 500
  }
  
  mean_init <- mean(working_df$counts)
  var_init <- var(working_df$counts)
  beta_init <- mean_init / var_init
  k_init <- mean_init * beta_init
  
  params_all <- getParameterCombinations(
    beta_range = seq(beta_init, max_range_multiple * beta_init, length.out = length_out),
    k_range = seq(k_init, max_range_multiple * k_init, length.out = length_out),
    lambda_range = lambda_range
  )
  
  print(paste0("########## Testing  ", nrow(params_all)," parameters #################"))
  
  # set max and min vals for weight at far limit of distribution
  if (!fix_weight) {
    weight_value <- 0.99 * max(working_df$counts)
  }
  
  #run CVM to find cost at each set of parameters
  cvm_stat_list <- plyr::ddply(.data = params_all, 
                               .variables = "params_num", 
                               .fun = getCVMDistance, 
                               "working_df" = working_df, 
                               "weight_value" = weight_value)
  
  # make a heat map of cvm_stat_list
  cvm_stat_list <- dplyr::inner_join(x = cvm_stat_list, y = params_all, by = "params_num")
  if (plot_data) {
    print(pbsR:::plotCVMHeatmap(cvm_stat_list))
  }
  
  #print("########## Distribution Parameters tested #################")
  #print(cvm_stat_list)
  
  # retrieve parameters at minimum cost
  optim_params <- params_all[which.min(cvm_stat_list$omega2), ]
  #optim_params$Name <- working_df$Name[1]
  return(optim_params)
}

#' @export
#' @importFrom magrittr %>%
getDistributionParametersWithOptim <- function(working_df, lambda_range = seq(0.6, 0.9, by = 0.1), length_out = 3,
                                               fix_weight = TRUE, weight_value = 500, plot_data = FALSE, plot_terra = FALSE,
                                               bin_width = 0.05, max_dens = 10) {
  counts_col <- grepl(pattern = "count", x = names(working_df), ignore.case = TRUE)
  if (sum(counts_col) != 1) {
    stop('Working_df must have exactly one column with a name similar to "counts".')
  }
  names(working_df)[counts_col] <- "counts"
  working_df <- working_df %>% dplyr::filter(counts > 0)
  params_init <- getDistributionParameters(working_df = working_df, lambda_range = lambda_range, length_out = length_out)
  # set max and min vals for weight at far limit of distribution
  beta_init <- params_init$beta
  k_init <- params_init$k
  lambda_init <- params_init$lambda
  if (!fix_weight) {
    weight_value <- 0.99 * max(working_df$counts)
  }
  # use optim to calculate parameters
  optim_params <- stats::optim(
    par = c(beta_init, k_init, lambda_init),
    fn = function(par) {
      params_df <- data.frame(beta = par[1], k = par[2], lambda = par[3])
      getCVMDistance(
        params_df = params_df,
        weight_value = weight_value, 
        working_df = working_df
      )
    },
    method = "L-BFGS-B", lower = c(1e-5, 1e-5, 1e-5), upper = c(Inf, Inf, 0.99)
  )
  # calculate lambda based on  p values
  lambda <- 2 * mean(stats::pgamma(
    q = working_df$counts, shape = optim_params$par[2],
    rate = optim_params$par[1], lower.tail = FALSE
  ) > 0.5)

  # get fit quality metric
  params_df <- data.frame(beta = optim_params$par[1], k = optim_params$par[2], lambda = lambda)
  # fit_quality <- findAreasWithPlot(working_df = working_df, params_df = params_df, xlim = max_dens, printGraph = TRUE)
  # params_df$resid_area <- fit_quality$RH
  return(params_df)
}
