library(pheatmap)
library(wesanderson)
library(dendextend)
library(psych)
args = commandArgs(trailingOnly=TRUE)
RNA_expData <- read.table(args[1], header = TRUE, row.names = 1)
V1 <- strsplit(args[2], ",")[[1]]
sampleStage <- strsplit(args[3], ",")[[1]]
V2 <- sampleStage
sampInfo<-data.frame(SampleID=V1,Group =V2)

sampInfo[,2] <- factor(sampInfo[,2],levels= sampleStage)
sampInfo <- sampInfo[order(sampInfo[,2]),]
matchIndexes <- match(sampInfo[,1],colnames(RNA_expData))
RNA_expData <- RNA_expData[,matchIndexes]
colnames(RNA_expData) <- as.character(sampInfo[,2])
RNA_expData <- RNA_expData[which(rowSums(RNA_expData) > 0),]

p <- as.dist(cor2dist(cor(RNA_expData)))
hc <- hclust(p)
p.dend <- reorder(as.dendrogram(hc),wts=order(match(sampleStage,colnames(RNA_expData))))
col_cluster <- as.hclust(p.dend)
pdf(file=args[4],width=10,height=15)
p <- pheatmap(RNA_expData,show_rownames=FALSE,scale="row",cluster_rows=FALSE,cluster_cols=col_cluster,fontsize=10,color = colorRampPalette(c("#6B70B0","white","red"))(100))
plot(p$tree_col)
dev.off()