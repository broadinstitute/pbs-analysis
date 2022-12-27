#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org/"), Ncpus = 8)

install.packages("BiocManager")
BiocManager::install("rtracklayer")