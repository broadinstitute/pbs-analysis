## code to prepare `hg19_5000_map_75` dataset goes here

hg19_5000_map_75 = read.csv("references/hg19_5000_map_75.bedgraph", sep = "\t", header=F)
colnames(hg19_5000_map_75) = c("chr","start","end","mappability_score")
hg19_5000_map_75$start = hg19_5000_map_75$start + 1
usethis::use_data(hg19_5000_map_75, overwrite = TRUE)
