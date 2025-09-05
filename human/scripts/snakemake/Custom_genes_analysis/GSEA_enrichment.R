library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genekitr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

file <- args[1]
geneIDs <- readLines(file)
expData <- read.csv(args[2],row.names = 1)
exprIDs <- rownames(expData)
exprIDs <- unlist(lapply(strsplit(exprIDs,"\\|"),function(x) x[[1]]))
exprIDs <- gsub("\\.\\d+","",exprIDs)
rownames(expData) <- exprIDs
matchIndexes <- match(geneIDs,rownames(expData))
expData <- expData[matchIndexes,]
expData <- expData[, -1]
expData <- na.omit(expData)
entrezIDs <- transId(
   id = rownames(expData),
   transTo = "entrez", org = "human", keepNA = FALSE
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
   OrgDb = org.Hs.eg.db,
   keyType = "ENTREZID",
   pvalueCutoff = 1,
   pAdjustMethod = "fdr",
   scoreType = "std",
)
gsea_result = res@result
write.csv(gsea_result,args[3])

res1 <- gseKEGG(
   geneList = geneList,
   organism = "hsa",  # Human KEGG code
   keyType = "ncbi-geneid",  # If using EntrezID, otherwise adjust as needed
   pvalueCutoff = 1,
   pAdjustMethod = "fdr",
   scoreType = "std"
)
gsea_result1 = res1@result
write.csv(gsea_result1,args[4])