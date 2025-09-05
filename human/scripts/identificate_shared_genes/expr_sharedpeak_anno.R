library(genekitr)
library(org.Hs.eg.db)
args <- commandArgs(trailingOnly = TRUE)
DOWNpeak <- read.csv(args[1],header = T)
DOWNENTREZIDs <- unique(DOWNpeak$geneId)
peakIDs2 <- transId(
   id = DOWNENTREZIDs,
   transTo = "ensembl", org = "human", keepNA = FALSE
)

expData <- read.table(file=args[2],sep=",",header=TRUE,row.names=1,stringsAsFactors=FALSE)
exprIDs2 <- rownames(expData)[which(expData$log2FoldChange < as.numeric(args[3]) & expData$padj < as.numeric(args[4]))]
exprIDs2 <- unlist(lapply(strsplit(exprIDs2,"\\|"),function(x) x[[1]]))
exprIDs2 <- gsub("\\.\\d+","",exprIDs2)
sharedDOWN_IDs <- intersect(peakIDs2[,2],exprIDs2)
peakIDs2_filtered <- peakIDs2[peakIDs2[,2] %in% sharedDOWN_IDs, ]
matched_rows <- match(peakIDs2_filtered[,1], DOWNpeak$geneId)

if (any(!is.na(matched_rows))) {
   peakIDs2_filtered$V4 <- DOWNpeak$V4[matched_rows[!is.na(matched_rows)]]
}
peak <- read.table(args[5], header = TRUE)
peakIDs2_filtered <- cbind(peakIDs2_filtered, peak[match(peakIDs2_filtered[,3], peak$name), c("chr", "start", "end")])
peakIDs2_filtered_sorted <- peakIDs2_filtered[, c(2, 4, 5, 6)]
peakIDs2_filtered_sorted <- peakIDs2_filtered_sorted[, c(2, 3,4, 1)]
write.table(peakIDs2_filtered_sorted,args[6], sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)

peakIDs1_filtered_sorted <- read.table(args[7], header = F)
names(peakIDs2_filtered_sorted) <- names(peakIDs1_filtered_sorted)
Sharedgene_bedfile <- rbind(peakIDs1_filtered_sorted,peakIDs2_filtered_sorted)
write.table(Sharedgene_bedfile,args[8], sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
