library(GenomicFeatures)
library(BSgenome.Taestivum.IWGSC.v1)
library(genekitr)

args <- commandArgs(trailingOnly = TRUE)
txdb_file <- args[3]
txdb <- loadDb(txdb_file)

Shared_UPgene_Inpromote <- read.table(args[1], header=FALSE, stringsAsFactors=FALSE)
GeneID <- Shared_UPgene_Inpromote$V1 
geneIDs <- transId(
   id = GeneID,
   transTo = "ensembl", org = "Taestivum", keepNA = FALSE
)
geneIDs <- na.omit(geneIDs)
ensemblIDs <- geneIDs$ensembl
tx_by_gene <- transcriptsBy(txdb, by = "gene")
valid_ensemblIDs <- ensemblIDs[ensemblIDs %in% names(tx_by_gene)]
transcriptCoordsByGene.GRangesList <- tx_by_gene[valid_ensemblIDs]

filtered_tx <- keepStandardChromosomes(transcriptCoordsByGene.GRangesList, pruning.mode="coarse")

filtered_tx <- endoapply(filtered_tx, function(gr) gr[start(gr) > 2000])
filtered_tx <- filtered_tx[elementNROWS(filtered_tx) > 0]

promoter.seqs <- getPromoterSeq(filtered_tx, Taestivum, upstream=2000, downstream=20)

promoter.seqs <- unlist(promoter.seqs)
writeXStringSet(promoter.seqs, args[2])