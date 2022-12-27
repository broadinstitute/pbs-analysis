#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org/"), Ncpus = 8)

# non bioconductor packages
install.packages("dplyr")
install.packages("optparse")
install.packages("plyr")
install.packages("magrittr")
install.packages("ggplot2")
install.packages("dip.test")
# Bioconductor packages
install.packages("BiocManager")
BiocManager::install("CNAnorm")
BiocManager::install("Rsubread")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("biovizBase")
