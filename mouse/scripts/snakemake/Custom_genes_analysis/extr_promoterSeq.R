library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(genekitr)
args <- commandArgs(trailingOnly = TRUE)
gene <- read.table(args[1],header=FALSE,stringsAsFactors=FALSE)
EnsembleId <- gene$V1 
geneIDs <- transId(
   id = EnsembleId,
   transTo = "entrez", org = "mouse", keepNA = FALSE
)
geneIDs <- na.omit(geneIDs)
entrezIDs <- geneIDs$entrezid
valid_entrezIDs <- entrezIDs[entrezIDs %in% names(transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene"))]
transcriptCoordsByGene.GRangesList <-
   transcriptsBy (TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")[valid_entrezIDs]

promoter.seqs <- getPromoterSeq(transcriptCoordsByGene.GRangesList,
                                     Mmusculus, upstream=2000, downstream=20)

promoter.seqs <- unlist(promoter.seqs)
writeXStringSet(promoter.seqs, args[2])