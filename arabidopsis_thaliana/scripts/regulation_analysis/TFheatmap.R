library(ComplexHeatmap)
library(circlize)
args <- commandArgs(trailingOnly = TRUE)
valscaling <- function(x){
   x <- (x/sum(x))
}

zscore <- function(x){
   x <- (x-mean(x))/sd(x)
}
Enriched_tfs <- read.table(file=args[1],header=FALSE,sep="\t",stringsAsFactors=FALSE)
Enriched_tfs <- Enriched_tfs[!duplicated(Enriched_tfs[,2]), ]
expData <- read.csv(args[2],row.names = 1)
exprIDs <- rownames(expData)
exprIDs <- unlist(lapply(strsplit(exprIDs,"\\|"),function(x) x[[1]]))
exprIDs <- gsub("\\.\\d+","",exprIDs)
rownames(expData) <- exprIDs
matchIndexes <- match(Enriched_tfs[,1],rownames(expData))
matchIndexes <- na.omit(matchIndexes)
matchIndexes <- unique(matchIndexes)
expData <- expData[matchIndexes,]
expData <- expData[!is.na(matchIndexes), ]
expData <- expData[, -1]
expData <- na.omit(expData)

UpTFs <- rownames(expData)[which(expData$log2FoldChange > as.numeric(args[3]) & expData$padj < as.numeric(args[4]))]   
DownTFs <- rownames(expData)[which(expData$log2FoldChange < as.numeric(args[5]) & expData$padj < as.numeric(args[6]))]
matchIndexes <- match(c(UpTFs,DownTFs),rownames(expData))
rownames(expData) <-  Enriched_tfs[match(rownames(expData),Enriched_tfs[,1]),2]
expData <- expData[, 7:10]

expData <- apply(expData,2,valscaling)
expData <- t(apply(expData,1,zscore))
expData <- na.omit(expData)
f1 = colorRamp2(seq(min(expData), max(expData), length = 4), c("#CCFFFF","#99CCFF","#FF9966", "#CC3333"))
pdf(file=args[7],width=3.5,height=15)
ha = rowAnnotation(foo = anno_mark(at = matchIndexes, labels = rownames(expData)[matchIndexes]))
Heatmap(expData, name = "z-score", cluster_columns = FALSE, show_row_names = FALSE, right_annotation = ha, col = f1, row_names_gp = gpar(fontsize = 4),row_km = 2)
dev.off()
