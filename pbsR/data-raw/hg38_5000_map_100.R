## code to prepare `hg38_5000_map_100` dataset goes here

hg38_5000_map_100 = read.csv("references/hg38_5000_map_100.bedgraph", sep = "\t", header=F)
colnames(hg38_5000_map_100) = c("chr","start","end","mappability_score")
hg38_5000_map_100$start = hg38_5000_map_100$start + 1
usethis::use_data(hg38_5000_map_100, overwrite = TRUE)
