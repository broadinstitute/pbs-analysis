#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(ggplot2))

# This source file contains functions used in different versions of the Fit-Quality Metric for evaluating the fit of a gamma distribution to
# the density of read counts for histone ChiP-seq data.
#
# Approaches used:
#   1) Looking at the difference between the CDF functions and seeing if there are patterns in the curves that could indicate a good fit froma bad one
#   2) Calculating the CVM and plotting extremeties to see whether there are connections between good fits and bad ones
#         - while this did not work, it did show that very very high CVMs tend to relate to quality of read depth
#   3)*** Find and sum incremental differences before and after the median based on qualities we expect in bad fits:
#         - LH fit: summing the absolute value difference between fitted and empirical until median
#         - RH fit: summing the difference between fitted and empirical after median IF fitted > empirical (overestimating amount of signal)
#
# Approach #3 was shown to be the best approach, and the functions below are used to calculate the fits.
# Code written by: Sreshtaa Rajesh

#Helper function that returns a list including x and y points that define the density of the Counts column of a working dataframe
getDensityPts <- function(working_df, xlim = 10) {
  densityPts <- density(working_df$counts, n = 1024, to = xlim)
  return(densityPts)
}

#Helper function that returns a list of points to make a gamma curve based on a given vector of x values and the gamma parameters
getFittedPts <- function(xValues, params_df) {
  fittedPts <- dgamma(x = xValues, shape = params_df$k, rate = params_df$beta)
  fittedPts <- as.data.frame(fittedPts)
  names(fittedPts) <- c("yFitted")
  return(fittedPts)
}

#Generates a one-row dataframe with the Left hand (LH) and Right hand (RH) areas. LH area is tha absolute value of the difference between
#the fitted and empirical curves up until the median. The RH area sums all the differences where the fitted curve is greater than the empirical
#curve after the median. NOTE: THE RH AREA IS A BETTER INDICATOR OF FIT
findAreasWithPlot <- function(working_df, params_df, xlim = 10, filetitle = NULL, printGraph = FALSE) { #switching printGraph to true will generate a plot similar to those from makePlot function
  densityPts <- getDensityPts(working_df = working_df, xlim = xlim)
  fittedPts <- getFittedPts(xValues = densityPts$x, params_df = params_df)
  lambda <- params_df$lambda
  difference_df <- data.frame(xValues = densityPts$x, yEmpirical = densityPts$y, yFitted = lambda*fittedPts)
  areas <- data.frame(LH = double(1), RH = double(1))

  #Right-Hand Area
  filtered_difference <- difference_df %>% dplyr::filter(xValues > median(working_df$counts)) %>% dplyr::filter(yEmpirical < yFitted)
  filtered_difference$difference <- (filtered_difference$yFitted - filtered_difference$yEmpirical)
  rhArea <- sum(filtered_difference$difference)

  #Left-Hand Area
  filtered_difference <- difference_df %>% dplyr::filter(xValues < median(working_df$counts))
  filtered_difference$difference <- abs(filtered_difference$yFitted - filtered_difference$yEmpirical)
  lhArea <- sum(filtered_difference$difference)

  areas[1,] <- c(lhArea,rhArea)

  #Plotting (only done if printGraph is specified to be true)
  if(printGraph) {
    diffPlot <- ggplot(data = difference_df, aes(x = xValues, y = yEmpirical, col = "Empirical")) + geom_line() + geom_line(data = difference_df, aes(y = yFitted, col = "Fitted"))
    diffPlot <- diffPlot + geom_vline(xintercept = median(working_df$counts), linetype = "dotted") + labs(x = "Counts", y = "Density")
    if(!is.null(title)) {
      diffPlot <- diffPlot + ggtitle(paste(gsub(pattern = "Alignment Post Processing", replacement = "APP", filetitle), "Density Curves"))
    }
    print(diffPlot)
  }
  return(areas)
}

#Takes one row of the lookup table and generates a density plot with any additions that are wanted
makePlot <- function(filepath, omegaValue = NA, filename = NULL, celltype = NULL) {
  cleaned_df <- convertToCleanedDfExistingWig(filepath)
  params <- getExistingParams(filepath)

  #**ONLY USE IF YOU WANT THE CVM VALUE DISPLAYED**
  #if(is.na(omegaValue)) {
  #  omegaValue <- findOmega(filename, working_df = cleaned_df, params_df = params)
  #}

  gammaValues <- generateGammaValues(df = cleaned_df, parameters = params)
  densityPlot <- plotDoubleDensity(filetitle = filename, cleaned_df = cleaned_df, gammaValues = gammaValues)

  #**USE WHICHEVER ONE YOU WANT DEPENDING ON WHAT DATA YOU WANT DISPLAYED**
  #+ annotate("text", x = 8, y = 0.5, label = paste("CVM:", round(omegaValue, 2)))
  #+ annotate("text", x = 8, y = 0.5, label = paste("Cell Type:\n", celltype))

  #UNCOMMENT WHEN WRITING TO A PDF
  print(densityPlot + annotate("text", x = 8, y = 0.5, label = paste("Cell Type:\n", celltype)))
}

#Wrapper for the makePlot function; makes plot from a specified row from a specified lookup table
#Uses values from the row of the lookup dataframe to fill in for the necessary values required by makePlot
#Found it to be very useful in general if you want to quickly generate plots
makePlotWIndex <- function(index, lookupdf) {
  filepath = lookupdf$Binned.Filename[index]
  filename = lookupdf$Name[index]
  celltype = lookupdf$Cell.Types[index]
  makePlot(filepath = filepath, filename = filename, celltype = celltype)
}


