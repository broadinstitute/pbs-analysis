## code to prepare `hg19_5000_map_100` dataset goes here

hg19_5000_map_100 = read.csv("references/hg19_5000_map_100.bedgraph", sep = "\t", header=F)
colnames(hg19_5000_map_100) = c("chr","start","end","mappability_score")
hg19_5000_map_100$start = hg19_5000_map_100$start + 1
usethis::use_data(hg19_5000_map_100, overwrite = TRUE)
