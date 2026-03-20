library(GenomicFeatures)
library(org.Taestivum.eg.db)
library(genekitr)
args <- commandArgs(trailingOnly = TRUE)
txdb_file <- args[3]
txdb <- loadDb(txdb_file)
gene_list <- read.table(args[1],header = F)
entrezIDs <- transId(
   id = gene_list$V1,
   transTo = "ensembl", org = "Taestivum", keepNA = FALSE
)
gene <- entrezIDs[,2]
write.table(gene,file = args[2],row.names = F,quote = F,col.names = F)