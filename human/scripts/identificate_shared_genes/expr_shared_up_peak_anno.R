library(genekitr)
library(org.Hs.eg.db)
args <- commandArgs(trailingOnly = TRUE)
UPpeak <- read.csv(args[1],header = T)
UPENTREZIDs <- unique(UPpeak$geneId)
peakIDs1 <- transId(
   id = UPENTREZIDs,
   transTo = "ensembl", org = "human", keepNA = FALSE
)

expData <- read.table(file=args[2],sep=",",header=TRUE,row.names=1,stringsAsFactors=FALSE)
exprIDs1 <- rownames(expData)[which(expData$log2FoldChange > as.numeric(args[3]) & expData$padj < as.numeric(args[4]))]
exprIDs1 <- unlist(lapply(strsplit(exprIDs1,"\\|"),function(x) x[[1]]))
exprIDs1 <- gsub("\\.\\d+","",exprIDs1)
sharedUp_IDs <- intersect(peakIDs1[,2],exprIDs1)
peakIDs1_filtered <- peakIDs1[peakIDs1[,2] %in% sharedUp_IDs, ]
matched_rows <- match(peakIDs1_filtered[,1], UPpeak$geneId)

if (any(!is.na(matched_rows))) {
   peakIDs1_filtered$V4 <- UPpeak$V4[matched_rows[!is.na(matched_rows)]]
}
peak <- read.table(args[5], sep = "\t", header = TRUE)
peakIDs1_filtered <- cbind(peakIDs1_filtered, peak[match(peakIDs1_filtered[,3], peak$name), c("chr", "start", "end")])
peakIDs1_filtered_sorted <- peakIDs1_filtered[, c(2, 4, 5, 6)]
peakIDs1_filtered_sorted <- peakIDs1_filtered_sorted[, c(2, 3,4, 1)]

write.table(peakIDs1_filtered_sorted,args[6], sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
