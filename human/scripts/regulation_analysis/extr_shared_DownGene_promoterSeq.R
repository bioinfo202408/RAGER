library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genekitr)
#DOwn
args <- commandArgs(trailingOnly = TRUE)
Shared_DOWNgene_Inpromote <- read.table(args[1],header=FALSE,stringsAsFactors=FALSE)
EnsembleId_down <- Shared_DOWNgene_Inpromote$V1 
geneIDs_down <- transId(
   id = EnsembleId_down,
   transTo = "entrez", org = "human", keepNA = FALSE
)
geneIDs_down <- na.omit(geneIDs_down)
entrezIDs_down <- geneIDs_down$entrezid
valid_entrezIDs_down <- entrezIDs_down[entrezIDs_down %in% names(transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))]
transcriptCoordsByGene.GRangesList_down <-
   transcriptsBy (TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")[valid_entrezIDs_down]

promoter.seqs_down <- getPromoterSeq(transcriptCoordsByGene.GRangesList_down,
                                     Hsapiens, upstream=2000, downstream=20)

promoter.seqs_down <- unlist(promoter.seqs_down)
writeXStringSet(promoter.seqs_down, args[2])