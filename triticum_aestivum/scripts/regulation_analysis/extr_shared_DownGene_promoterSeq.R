library(GenomicFeatures)
library(BSgenome.Taestivum.IWGSC.v1)
library(genekitr)

args <- commandArgs(trailingOnly = TRUE)
txdb_file <- args[3]
txdb <- loadDb(txdb_file)

Shared_DOWNgene_Inpromote <- read.table(args[1], header = FALSE, stringsAsFactors = FALSE)
GeneID_down <- Shared_DOWNgene_Inpromote$V1

geneIDs_down <- transId(
   id = GeneID_down,
   transTo = "ensembl", org = "Taestivum", keepNA = FALSE
)
geneIDs_down <- na.omit(geneIDs_down)

ensemblIDs_down <- geneIDs_down$ensembl
valid_ensemblIDs_down <- ensemblIDs_down[ensemblIDs_down %in% names(transcriptsBy(txdb, by = "gene"))]

transcriptCoordsByGene.GRangesList_down <-
   transcriptsBy(txdb, by = "gene")[valid_ensemblIDs_down]

promoter.seqs_down <- getPromoterSeq(
  transcriptCoordsByGene.GRangesList_down,
  Taestivum,
  upstream = 2000,
  downstream = 20
)

promoter.seqs_down <- unlist(promoter.seqs_down)
writeXStringSet(promoter.seqs_down, args[2])