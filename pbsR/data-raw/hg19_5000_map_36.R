## code to prepare `hg19_5000_map_36` dataset goes here

hg19_5000_map_36 = read.csv("references/hg19_5000_map_36.bedgraph", sep = "\t", header=F)
colnames(hg19_5000_map_36) = c("chr","start","end","mappability_score")
hg19_5000_map_36$start = hg19_5000_map_36$start + 1
usethis::use_data(hg19_5000_map_36, overwrite = TRUE)
