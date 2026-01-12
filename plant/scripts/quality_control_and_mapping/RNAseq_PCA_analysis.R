library(ggplot2)
library(scatterplot3d)

args <- commandArgs(trailingOnly = TRUE)

expData <- read.table(file = args[1], sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

V1 <- strsplit(args[2], ",")[[1]]  
sampleStage <- strsplit(args[3], ",")[[1]]  
V2 <- sampleStage
sampInfo <- data.frame(SampleID = V1, Group = V2, stringsAsFactors = FALSE)

group_names <- strsplit(args[4], ",")[[1]]  
sample_counts <- as.numeric(strsplit(args[5], ",")[[1]])  

sampInfo[, 2] <- factor(sampInfo[, 2], levels = sampleStage)
sampInfo <- sampInfo[order(sampInfo[, 2]), ]
matchIndexes <- match(sampInfo[, 1], colnames(expData))
expData <- expData[, matchIndexes]
colnames(expData) <- as.character(sampInfo[, 2])

stageNames <- rep(group_names, sample_counts)

PCAFitClass <- prcomp(t(expData))
max3PCs <- PCAFitClass$x[, 1:3]

group_factor <- factor(stageNames, levels = group_names)
colors <- rainbow(length(group_names))[group_factor]
shapes <- seq(1, length(group_names))[group_factor]

plot_PCA <- function(x_pc, y_pc, xlabel, ylabel, outfile) {
   x <- (x_pc - min(x_pc)) / (max(x_pc) - min(x_pc))
   y <- (y_pc - min(y_pc)) / (max(y_pc) - min(y_pc))
   plotdata <- data.frame(x = x, y = y, pointColor = stageNames)
   
   p <- ggplot(plotdata, aes(
      x = x, y = y,
      colour = factor(pointColor, levels = group_names),
      shape = factor(pointColor, levels = group_names)
   )) +
      geom_point(size = 5) +
      scale_shape_manual(values = seq(1, length(group_names))) +
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

plot_PCA(max3PCs[, 1], max3PCs[, 2], "PC1", "PC2", args[6])
plot_PCA(max3PCs[, 1], max3PCs[, 3], "PC1", "PC3", args[7])
plot_PCA(max3PCs[, 2], max3PCs[, 3], "PC2", "PC3", args[8])

x <- (max3PCs[, 1] - min(max3PCs[, 1])) / (max(max3PCs[, 1]) - min(max3PCs[, 1]))
y <- (max3PCs[, 2] - min(max3PCs[, 2])) / (max(max3PCs[, 2]) - min(max3PCs[, 2]))
z <- (max3PCs[, 3] - min(max3PCs[, 3])) / (max(max3PCs[, 3]) - min(max3PCs[, 3]))

pdf(file = args[9], width = 18, height = 9)
par(xpd = TRUE, mgp = c(2, 0.8, 0), oma = c(0.2, 0.2, 0.2, 0.2), lwd = 3)
scatterplot3d(
   x, y, z,
   color = colors,
   pch = shapes,
   xlab = "PC1", ylab = "PC2", zlab = "PC3",
   mar = c(4, 4, 4, 28),
   cex.axis = 2, cex.lab = 3, cex.symbols = 2.5
)
legend("topright",
       legend = group_names,
       inset = c(-0.47, -0.06),
       pch = seq(1, length(group_names)),
       col = rainbow(length(group_names)),
       cex = 2, y.intersp = 1.0
)
dev.off()
