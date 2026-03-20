args <- commandArgs(trailingOnly = TRUE)
gene_anno <- read.table(args[1],header=FALSE,stringsAsFactors=FALSE)
gene_anno$V4 <- sub("\\..*$", "", gene_anno$V4)
expData <- read.table(file=args[2],sep=",",header=TRUE,row.names=1,stringsAsFactors=FALSE)
expData_UPgene <- expData[which((expData$log2FoldChange) > as.numeric(args[3]) & expData$padj < as.numeric(args[4])),]
exprIDs <- rownames(expData_UPgene)
exprIDs <- unlist(lapply(strsplit(exprIDs,"\\|"),function(x) x[[1]]))
exprIDs <- gsub("\\.\\d+","",exprIDs)

gene_anno_UPgene <- gene_anno[gene_anno$V4 %in% exprIDs, ]
write.table(gene_anno_UPgene, file = args[5], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

expData_DOWNgene <- expData[which((expData$log2FoldChange) < as.numeric(args[6]) & expData$padj < as.numeric(args[7])),]
exprIDs1 <- rownames(expData_DOWNgene)
exprIDs1 <- unlist(lapply(strsplit(exprIDs1,"\\|"),function(x) x[[1]]))
exprIDs1 <- gsub("\\.\\d+","",exprIDs1)

gene_anno_DOWNgene <- gene_anno[gene_anno$V4 %in% exprIDs1, ]
write.table(gene_anno_DOWNgene, file = args[8], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)