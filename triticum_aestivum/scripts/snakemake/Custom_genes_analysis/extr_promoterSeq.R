library(GenomicFeatures)
library(BSgenome.Taestivum.IWGSC.v1)
library(org.Taestivum.eg.db)
library(genekitr)
args <- commandArgs(trailingOnly = TRUE)
txdb_file <- args[3]
txdb <- loadDb(txdb_file)
gene <- read.table(args[1],header=FALSE,stringsAsFactors=FALSE)
EnsembleId <- gene$V1 
geneIDs <- transId(
   id = EnsembleId,
   transTo = "entrez", org = "Taestivum", keepNA = FALSE
)
geneIDs <- na.omit(geneIDs)
entrezIDs <- geneIDs$entrezid
valid_entrezIDs <- entrezIDs[entrezIDs %in% names(transcriptsBy(txdb, by = "gene"))]
transcriptCoordsByGene.GRangesList <-
   transcriptsBy (txdb, by = "gene")[valid_entrezIDs]

promoter.seqs <- getPromoterSeq(transcriptCoordsByGene.GRangesList,
                                     Taestivum, upstream=2000, downstream=20)

promoter.seqs <- unlist(promoter.seqs)
writeXStringSet(promoter.seqs, args[2])