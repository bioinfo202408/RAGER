library(ggplot2)
library(scatterplot3d)

args <- commandArgs(trailingOnly = TRUE)


input_file <- args[1]  
col_ranges <- strsplit(args[2], ",")[[1]]  
colnames_list <- strsplit(args[3], ",")[[1]]  
sample_names <- strsplit(args[4], ",")[[1]]  
sample_counts <- as.numeric(strsplit(args[5], ",")[[1]])  
out_PC1_PC2 <- args[6]
out_PC1_PC3 <- args[7]
out_PC2_PC3 <- args[8]
out_3D <- args[9]


atacData_raw <- read.table(file = input_file, header = FALSE, stringsAsFactors = FALSE, skip = 3)


sample_means <- list()
for (range in col_ranges) {
   idx <- as.numeric(eval(parse(text = range)))
   idx <- idx[idx <= ncol(atacData_raw)]  
   if (length(idx) == 0) {
      stop(paste("Column range", range, "Invalid in the data"))
   }
   sample_means[[length(sample_means) + 1]] <- rowMeans(atacData_raw[, idx, drop = FALSE], na.rm = TRUE)
}


atacData <- do.call(cbind, sample_means)
colnames(atacData) <- colnames_list
stageNames <- rep(sample_names, sample_counts)


atacData <- atacData[complete.cases(atacData) & !apply(atacData, 1, function(x) any(is.infinite(x))), ]
atacData[is.na(atacData)] <- 0


PCAFitClass <- prcomp(t(atacData))
max3PCs <- PCAFitClass$x[, 1:3]


plot_PCA_2D <- function(x_pc, y_pc, xlabel, ylabel, outfile) {
   x <- (x_pc - min(x_pc)) / (max(x_pc) - min(x_pc))
   y <- (y_pc - min(y_pc)) / (max(y_pc) - min(y_pc))
   plotdata <- data.frame(x = x, y = y, pointColor = stageNames)
   
   p <- ggplot(plotdata, aes(
      x = x, y = y,
      colour = factor(pointColor, levels = sample_names),
      shape = factor(pointColor, levels = sample_names)
   )) +
      geom_point(size = 5) +
      scale_shape_manual(values = seq(1, length(sample_names))) +
      xlab(xlabel) + ylab(ylabel) +
      theme(
         legend.title = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = "black"),
         axis.title.x = element_text(size = 14),
         axis.title.y = element_text(size = 14),
         axis.text.x = element_text(size = 11),
         axis.text.y = element_text(size = 11)
      )
   ggsave(outfile, plot = p, width = 8, height = 6)
}


plot_PCA_2D(max3PCs[, 1], max3PCs[, 2], "PC1", "PC2", out_PC1_PC2)
plot_PCA_2D(max3PCs[, 1], max3PCs[, 3], "PC1", "PC3", out_PC1_PC3)
plot_PCA_2D(max3PCs[, 2], max3PCs[, 3], "PC2", "PC3", out_PC2_PC3)


x <- (max3PCs[, 1] - min(max3PCs[, 1])) / (max(max3PCs[, 1]) - min(max3PCs[, 1]))
y <- (max3PCs[, 2] - min(max3PCs[, 2])) / (max(max3PCs[, 2]) - min(max3PCs[, 2]))
z <- (max3PCs[, 3] - min(max3PCs[, 3])) / (max(max3PCs[, 3]) - min(max3PCs[, 3]))

group_factor <- factor(stageNames, levels = sample_names)
colors <- rainbow(length(sample_names))[group_factor]
shapes <- seq(1, length(sample_names))[group_factor]

pdf(file = out_3D, width = 18, height = 9)
par(xpd = TRUE, mgp = c(2, 0.8, 0), oma = c(0.2, 0.2, 0.2, 0.2), lwd = 3)
scatterplot3d(
   x, y, z,
   color = colors,
   pch = shapes,
   xlab = "PC1", ylab = "PC2", zlab = "PC3",
   mar = c(4, 4, 4, 28),
   cex.axis = 2, cex.lab = 3, cex.symbols = 5
)
legend("topright", legend = sample_names,
       inset = c(-0.47, -0.06),
       pch = seq(1, length(sample_names)),
       col = rainbow(length(sample_names)),
       cex = 2, y.intersp = 1.0
)
dev.off()
