source("./map.enhancers.to.genes.R")
options(echo=T)
args <- commandArgs(trailingOnly = TRUE)
targets <- args[1]
interactions <- args[2]

annotated.enhancers <- map.enhancers.to.genes(targets, interactions)

write.table(annotated.enhancers, gsub("$", ".annotated.bed", targets), sep = "\t", quote = F, col.names = F, row.names = F)
