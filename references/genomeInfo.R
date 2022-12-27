genomes <- c('BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg18', 
             'BSgenome.Mmusculus.UCSC.mm10', 'BSgenome.Mmusculus.UCSC.mm9')

for (genome in genomes){
  name <- strsplit(genome, '\\.')[[1]][4]
  
  library(genome, character.only=T)
  
  genomeInfo <- seqinfo(eval(parse(text=genome)))
  
  save(genomeInfo, file=sprintf('./references/%s.RData', name))
}
