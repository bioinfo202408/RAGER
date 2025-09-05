library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genekitr)
args <- commandArgs(trailingOnly = TRUE)
gene <- read.table(args[1],header=FALSE,stringsAsFactors=FALSE)
EnsembleId <- gene$V1 
geneIDs <- transId(
   id = EnsembleId,
   transTo = "entrez", org = "human", keepNA = FALSE
)
geneIDs <- na.omit(geneIDs)
entrezIDs <- geneIDs$entrezid
valid_entrezIDs <- entrezIDs[entrezIDs %in% names(transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))]
transcriptCoordsByGene.GRangesList <-
   transcriptsBy (TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")[valid_entrezIDs]

promoter.seqs <- getPromoterSeq(transcriptCoordsByGene.GRangesList,
                                     Hsapiens, upstream=2000, downstream=20)

promoter.seqs <- unlist(promoter.seqs)
writeXStringSet(promoter.seqs, args[2])