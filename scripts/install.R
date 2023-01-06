#!/usr/bin/env Rscript

# Use all available CPU cores for the build
library(parallel)
options(Ncpus = detectCores())

# helper to install and load packages one by one
# to verify they were installed successfully;
# otherwise, R treats installation errors as warnings
setup <- function(installer, packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      installer(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# Install randomForest directly from url because installer can't find it otherwise
rf <- 'https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz'
install.packages(rf, repos=NULL, type="source") 
library(randomForest)

setup(install.packages, c(
  # for rescaling
  "data.table",
  "dplyr",
  "optparse",

  # for cnv rescaling
  "diptest",
  "ggplot2",
  "magrittr",
  "plyr",

  # for fitting
  "fBasics",
  "goftest",
  "reshape",
  "tidyr",

  # for Bioconductor packages
  "BiocManager"
))

setup(BiocManager::install, c(
  "rtracklayer",

  # for cnv rescaling
  #"biovizBase",
  #"BSgenome.Hsapiens.UCSC.hg19",
  #"BSgenome.Hsapiens.UCSC.hg38",
  "CNAnorm"
  #"Rsubread"
))
