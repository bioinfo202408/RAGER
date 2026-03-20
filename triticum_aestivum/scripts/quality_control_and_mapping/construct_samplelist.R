args <- commandArgs(trailingOnly = TRUE)
V1 <- strsplit(args[1], ",")[[1]]  
V2 <- strsplit(args[2], ",")[[1]]
sampInfo <- data.frame(SampleID = V1, Group = V2, stringsAsFactors = FALSE)
write.table(sampInfo, args[3], sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
