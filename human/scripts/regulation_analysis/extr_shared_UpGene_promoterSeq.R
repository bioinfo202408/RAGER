library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genekitr)
args <- commandArgs(trailingOnly = TRUE)
Shared_UPgene_Inpromote <- read.table(args[1], header=FALSE, stringsAsFactors=FALSE)
EnsembleId <- Shared_UPgene_Inpromote$V1 
geneIDs <- transId(
   id = EnsembleId,
   transTo = "entrez", org = "human", keepNA = FALSE
)
geneIDs <- na.omit(geneIDs)
entrezIDs <- geneIDs$entrezid
tx_by_gene <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
valid_entrezIDs <- entrezIDs[entrezIDs %in% names(tx_by_gene)]
transcriptCoordsByGene.GRangesList <- tx_by_gene[valid_entrezIDs]

filtered_tx <- keepStandardChromosomes(transcriptCoordsByGene.GRangesList, pruning.mode="coarse")

filtered_tx <- endoapply(filtered_tx, function(gr) gr[start(gr) > 2000])
filtered_tx <- filtered_tx[elementNROWS(filtered_tx) > 0]

promoter.seqs <- getPromoterSeq(filtered_tx, Hsapiens, upstream=2000, downstream=20)

promoter.seqs <- unlist(promoter.seqs)
writeXStringSet(promoter.seqs, args[2])


