args <- commandArgs(trailingOnly = TRUE)
enh_With1Mb <- read.table(args[1],header=FALSE,stringsAsFactors=FALSE)
unique_UPenh <- unique(enh_With1Mb$V1)
enh_anno <- read.table(args[2],header=FALSE,stringsAsFactors=FALSE)
UPenh_anno <- enh_anno[enh_anno$V4 %in% unique_UPenh, ]
write.table(UPenh_anno, file = args[3],  sep = "\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

DOWNenh_With1Mb <- read.table(args[4],header=FALSE,stringsAsFactors=FALSE)
unique_DOWNenh <- unique(DOWNenh_With1Mb$V1)
DOWNenh_anno <- enh_anno[enh_anno$V4 %in% unique_DOWNenh, ]
write.table(DOWNenh_anno, file = args[5], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
