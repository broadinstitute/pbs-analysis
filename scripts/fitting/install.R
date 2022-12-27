#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org/"), Ncpus = 8)

install.packages("dplyr")
install.packages("reshape")
install.packages("ggplot2")
install.packages("fBasics")
install.packages("goftest")
install.packages("tidyr")
suppressPackageStartupMessages(library(reshape))
