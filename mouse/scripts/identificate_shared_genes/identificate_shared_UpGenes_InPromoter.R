library(ggplot2)
library(genekitr)
library(ggVennDiagram)
library(org.Mm.eg.db)
args <- commandArgs(trailingOnly = TRUE)
UPpeak <- read.csv(args[1],header = T)

UPENTREZIDs <- unique(UPpeak$geneId)
peakIDs1 <- transId(
   id = UPENTREZIDs,
   transTo = "ensembl", org = "mouse", keepNA = FALSE
)

expData <- read.table(file=args[2],sep=",",header=TRUE,row.names=1,stringsAsFactors=FALSE)
exprIDs1 <- rownames(expData)[which(expData$log2FoldChange > as.numeric(args[3]) & expData$padj < as.numeric(args[4]))]
exprIDs1 <- unlist(lapply(strsplit(exprIDs1,"\\|"),function(x) x[[1]]))
exprIDs1 <- gsub("\\.\\d+","",exprIDs1)

plotData <- list(peak=peakIDs1[,2],gene=exprIDs1)
pdf(file=args[5],width=8,height=6)
ggVennDiagram(plotData, 
              edge_size = 0.5) + 
   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
   ggtitle("shared up genes") +
   theme(plot.title = element_text(hjust = 0.5))
dev.off()
sharedUp_IDs <- intersect(peakIDs1[,2],exprIDs1)
write.table(sharedUp_IDs,file=args[6],row.names=FALSE,col.names=FALSE,quote=FALSE)

