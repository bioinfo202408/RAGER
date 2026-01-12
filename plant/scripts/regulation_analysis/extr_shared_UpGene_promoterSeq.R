library(TxDb.Athaliana.BioMart.plantsmart28)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(genekitr)

args <- commandArgs(trailingOnly = TRUE)
Shared_UPgene_Inpromote <- read.table(args[1], header=FALSE, stringsAsFactors=FALSE)
GeneID <- Shared_UPgene_Inpromote$V1 
geneIDs <- transId(
   id = GeneID,
   transTo = "ensembl", org = "arabidopsis", keepNA = FALSE
)
geneIDs <- na.omit(geneIDs)
ensemblIDs <- geneIDs$ensembl
tx_by_gene <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
valid_ensemblIDs <- ensemblIDs[ensemblIDs %in% names(tx_by_gene)]
transcriptCoordsByGene.GRangesList <- tx_by_gene[valid_ensemblIDs]

filtered_tx <- keepStandardChromosomes(transcriptCoordsByGene.GRangesList, pruning.mode="coarse")

filtered_tx <- endoapply(filtered_tx, function(gr) gr[start(gr) > 2000])
filtered_tx <- filtered_tx[elementNROWS(filtered_tx) > 0]

promoter.seqs <- getPromoterSeq(filtered_tx, Mmusculus, upstream=2000, downstream=20)

promoter.seqs <- unlist(promoter.seqs)
writeXStringSet(promoter.seqs, args[2])


