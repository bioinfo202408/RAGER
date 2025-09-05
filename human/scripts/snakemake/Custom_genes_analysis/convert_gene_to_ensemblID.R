library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(genekitr)
args <- commandArgs(trailingOnly = TRUE)
gene_list <- read.table(args[1],header = F)
entrezIDs <- transId(
   id = gene_list$V1,
   transTo = "ensembl", org = "human", keepNA = FALSE
)
gene <- entrezIDs[,2]
write.table(gene,file = args[2],row.names = F,quote = F,col.names = F)

            