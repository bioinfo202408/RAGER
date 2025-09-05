library(clusterProfiler)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(genekitr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

file <- args[1]
geneIDs <- readLines(file)
geneIDs <- unique(geneIDs)
expData <- read.csv(args[2],row.names = 1)
exprIDs <- rownames(expData)
exprIDs <- unlist(lapply(strsplit(exprIDs,"\\|"),function(x) x[[1]]))
exprIDs <- gsub("\\.\\d+","",exprIDs)
rownames(expData) <- exprIDs
matchIndexes <- match(geneIDs,rownames(expData))
expData <- expData[matchIndexes,]
expData <- expData[, -1]
entrezIDs <- transId(
   id = rownames(expData),
   transTo = "entrez", org = "mouse", keepNA = FALSE
)
entrezIDs <- entrezIDs[!duplicated(entrezIDs$input_id), ]
expData$input_id <- rownames(expData)

mergedData <- merge(entrezIDs, expData, by = "input_id")
sortedData <- mergedData[order(mergedData$log2FoldChange, decreasing = TRUE), ]
geneList <- sortedData$log2FoldChange
names(geneList) <- sortedData$entrezid
res <- gseGO(
   geneList,
   ont = "ALL",
   OrgDb = org.Mm.eg.db,
   keyType = "ENTREZID",
   pvalueCutoff = 1,
   pAdjustMethod = "fdr",
   scoreType = "neg",
)
gsea_result = res@result
write.csv(gsea_result,args[3])

res1 <- gseKEGG(
   geneList = geneList,
   organism = "mmu",
   keyType = "ncbi-geneid",
   pvalueCutoff = 1,
   pAdjustMethod = "fdr",
   scoreType = "neg",
)
gsea_result1 = res1@result
write.csv(gsea_result1,args[4])
saveRDS(res, args[5])
saveRDS(res1, args[6])