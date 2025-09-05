library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)

database_all <- read.table(file = args[1], sep = ",", header = TRUE, check.names = FALSE)
expdata <- database_all[1:(nrow(database_all) - 1), 2:ncol(database_all)]
rownames(expdata) <- database_all[1:(nrow(database_all) - 1), 1]

group_info <- strsplit(args[2], ",")[[1]]      
group_levels <- strsplit(args[3], ",")[[1]]  
group_info <- factor(group_info, levels = group_levels)

coldata <- data.frame(row.names = colnames(expdata), group = group_info)

dds <- DESeqDataSetFromMatrix(countData = expdata,
                              colData = coldata,
                              design = ~ group)
dds <- DESeq(dds)
res <- results(dds)

res <- res[order(res$padj, na.last = TRUE), ]
resdata <- merge(as.data.frame(res),
                 as.data.frame(counts(dds, normalized = TRUE)),
                 by = "row.names", sort = FALSE)
rownames(resdata) <- resdata[, 1]
write.csv(resdata, file = args[4], row.names = TRUE)