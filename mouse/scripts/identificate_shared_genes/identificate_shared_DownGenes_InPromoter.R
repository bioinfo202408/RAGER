library(ggplot2)
library(genekitr)
library(ggVennDiagram)
library(org.Mm.eg.db)
args <- commandArgs(trailingOnly = TRUE)
Downpeak <- read.csv(args[1],header = T)

DOWNENTREZIDs <- unique(Downpeak$geneId)
peakIDs2 <- transId(
   id = DOWNENTREZIDs,
   transTo = "ensembl", org = "mouse", keepNA = FALSE
)

expData <- read.table(file=args[2],sep=",",header=TRUE,row.names=1,stringsAsFactors=FALSE)
exprIDs2 <- rownames(expData)[which(expData$log2FoldChange < as.numeric(args[3]) & expData$padj < as.numeric(args[4]))]
exprIDs2 <- unlist(lapply(strsplit(exprIDs2,"\\|"),function(x) x[[1]]))
exprIDs2 <- gsub("\\.\\d+","",exprIDs2)

plotData2 <- list(Peak=peakIDs2[,2],Gene=exprIDs2)
pdf(file=args[5],width=8,height=6)
ggVennDiagram(plotData2) + 
   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
   ggtitle("shared down genes") +
   theme(plot.title = element_text(hjust = 0.5))
dev.off()
sharedDown_IDs <- intersect(peakIDs2[,2],exprIDs2)
write.table(sharedDown_IDs,file=args[6],row.names=FALSE,col.names=FALSE,quote=FALSE)
