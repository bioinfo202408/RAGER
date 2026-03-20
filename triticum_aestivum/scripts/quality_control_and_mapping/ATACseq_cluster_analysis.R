library(pheatmap)
library(wesanderson)
library(dendextend)
library(psych)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_clustering <- args[2]
output_reorder <- args[3]
sample_names <- strsplit(args[4], ",")[[1]]
col_ranges <- strsplit(args[5], ",")[[1]]
order_vec <- as.integer(strsplit(args[6], ",")[[1]])

atacData <- read.table(file = input_file, header = FALSE, stringsAsFactors = FALSE, skip = 3)

sample_means <- list()
for (range in col_ranges) {
   idx <- as.numeric(eval(parse(text = range)))
   idx <- idx[idx <= ncol(atacData)]
   if (length(idx) == 0) {
      stop(paste("Column range", range, "Invalid in the data"))
   }
   sample_means[[length(sample_means) + 1]] <- rowMeans(atacData[, idx, drop = FALSE], na.rm = TRUE)
}

mat <- do.call(cbind, sample_means)

mat <- mat[complete.cases(mat) & !apply(mat, 1, function(x) any(is.infinite(x))), ]

atacData <- do.call(cbind, sample_means)
colnames(atacData) <- sample_names
atacData <- atacData[rowSums(atacData) > 0,]

pdf(file = output_clustering, width = 10, height = 15)
p <- pheatmap(atacData, show_rownames = FALSE, scale = "row",
              cluster_rows = FALSE, cluster_cols = TRUE,
              clustering_distance_cols = "correlation",
              fontsize = 10, color = colorRampPalette(c("#6B70B0","white","red"))(100))
plot(p$tree_col)
dev.off()

col_cluster <- p$tree_col
col_cluster$order <- order_vec

pdf(file = output_reorder, width = 10, height = 15)
p <- pheatmap(atacData, show_rownames = FALSE, scale = "row",
              cluster_rows = FALSE, cluster_cols = col_cluster,
              fontsize = 10,
              color = colorRampPalette(c("#6B70B0","white","red"))(100))
plot(p$tree_col)
dev.off()
