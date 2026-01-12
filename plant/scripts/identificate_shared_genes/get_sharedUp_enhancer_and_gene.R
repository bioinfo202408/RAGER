args <- commandArgs(trailingOnly = TRUE)
annoData <- read.table(file=args[1],header=FALSE)
sortedRegions <- read.table(file=args[2],header=FALSE)

matchedColumn <- vector("numeric", length = nrow(sortedRegions))

for (i in 1:nrow(sortedRegions)) {
   regionValue <- sortedRegions[i, 2]
   matchingRow <- which(annoData[, 2] == regionValue)
   if (length(matchingRow) > 0) {
      matchedColumn[i] <- annoData[matchingRow, 4]
   } else {
      matchedColumn[i] <- NA
   }
}
sortedRegions$matchedColumn <- matchedColumn

atacDataraw <- read.table(file = args[3], header = FALSE, stringsAsFactors = FALSE, skip = 3)
atacDataraw[which(is.na(atacDataraw),arr.ind=TRUE)] <- 0
col_ranges <- strsplit(args[4], ",")[[1]]
sample_means <- list()
for (range in col_ranges) {
   idx <- as.numeric(eval(parse(text = range)))
   idx <- idx[idx <= ncol(atacDataraw)]
   if (length(idx) == 0) {
      stop(paste("col_range", range, "Invalid in the data"))
   }
   sample_means[[length(sample_means) + 1]] <- rowMeans(atacDataraw[, idx, drop = FALSE], na.rm = TRUE)
}

atacData <- do.call(cbind, sample_means)

colnames(atacData) <- strsplit(args[5], ",")[[1]]
rownames(atacData) <- sortedRegions$matchedColumn
atacData <- as.data.frame(atacData)
atacData$ATACA <- rowMeans(atacData[, strsplit(args[6], ",")[[1]]])
atacData$ATACN <- rowMeans(atacData[, strsplit(args[7], ",")[[1]]])

pseudoCount <- 1e-6

ATACA_adj <- ifelse(atacData$ATACA == 0, pseudoCount, atacData$ATACA)
ATACN_adj <- ifelse(atacData$ATACN == 0, pseudoCount, atacData$ATACN)

atacData$log2FoldChange <- log2(ATACA_adj / ATACN_adj)

filtered_data <- atacData[atacData$log2FoldChange > 1, ]
enh_With1Mb <- read.table(args[8],header=FALSE,stringsAsFactors=FALSE)

filtered_row_names <- rownames(filtered_data)
write.table(filtered_row_names,file=args[9], row.names = FALSE, col.names = FALSE, quote = FALSE)
enh_With1Mb_filtered <- enh_With1Mb[enh_With1Mb$V1 %in% filtered_row_names, ]

unique_values <- unique(enh_With1Mb_filtered$V2)
write.table(unique_values,file=args[10], row.names = FALSE, col.names = FALSE, quote = FALSE)

