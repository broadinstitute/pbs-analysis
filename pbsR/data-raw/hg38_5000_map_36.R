## code to prepare `hg38_5000_map_36` dataset goes here

hg38_5000_map_36 = read.csv("references/hg38_5000_map_36.bedgraph", sep = "\t", header=F)
colnames(hg38_5000_map_36) = c("chr","start","end","mappability_score")
hg38_5000_map_36$start = hg38_5000_map_36$start + 1
usethis::use_data(hg38_5000_map_36, overwrite = TRUE)
