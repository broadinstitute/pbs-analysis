## code to prepare `mm10_5000_map_100` dataset goes here

mm10_5000_map_100 = read.csv("data-raw/mm10_5000_map_100.bedgraph.bed", sep = "\t", header=F)
colnames(mm10_5000_map_100) = c("chr","start","end","mappability_score")
mm10_5000_map_100$start = mm10_5000_map_100$start + 1
usethis::use_data(mm10_5000_map_100, overwrite = TRUE)
