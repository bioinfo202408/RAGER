library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
expData <- read.table(file=args[1],sep=",",header=TRUE,row.names=1,stringsAsFactors=FALSE)

expData <- expData[which(abs(expData$log2FoldChange) > as.numeric(args[2]) & expData$padj < as.numeric(args[3])),]
exprIDs <- rownames(expData)
exprIDs <- unlist(lapply(strsplit(exprIDs,"\\|"),function(x) x[[1]]))
exprIDs <- gsub("\\.\\d+","",exprIDs)
rownames(expData) <- exprIDs

annoData <- read.table(file=args[4],header=FALSE)
annoData[,4] <- gsub("\\.\\d+","",annoData[,4])

atacData_raw <- read.table(file = args[5], header = FALSE, stringsAsFactors = FALSE, skip = 3)
atacData_raw[which(is.na(atacData_raw),arr.ind=TRUE)] <- 0
col_ranges <- strsplit(args[6], ",")[[1]]
sample_means <- list()
for (range in col_ranges) {
   idx <- as.numeric(eval(parse(text = range)))
   idx <- idx[idx <= ncol(atacData_raw)]
   if (length(idx) == 0) {
      stop(paste("col_range", range, "Invalid in the data"))
   }
   sample_means[[length(sample_means) + 1]] <- rowMeans(atacData_raw[, idx, drop = FALSE], na.rm = TRUE)
}

atacData <- do.call(cbind, sample_means)

colnames(atacData) <- strsplit(args[7], ",")[[1]]
rownames(atacData) <- annoData[, 4]
atacData <- as.data.frame(atacData)
atacData$ATACAdata <- rowMeans(atacData[, strsplit(args[8], ",")[[1]]])
atacData$ATACNdata <- rowMeans(atacData[, strsplit(args[9], ",")[[1]]])
atacData$log2FoldChange <- log2(atacData$ATACAdata / atacData$ATACNdata)

UpGenes <- read.table(file=args[10], header=FALSE, stringsAsFactors=FALSE)
DownGenes <- read.table(file=args[11], header=FALSE, stringsAsFactors=FALSE)

UpGenesMatched <- intersect(UpGenes$V1, row.names(atacData))
DownGenesMatched <- intersect(DownGenes$V1, row.names(atacData))

UpATACLog2FC <- atacData[UpGenesMatched, "log2FoldChange"]
DownATACLog2FC <- atacData[DownGenesMatched, "log2FoldChange"]

UpExpLog2FC <- expData[UpGenesMatched, "log2FoldChange"]
DownExpLog2FC <- expData[DownGenesMatched, "log2FoldChange"]

plotData <- data.frame(
   ExpLog2FC = c(UpExpLog2FC, DownExpLog2FC),
   ATACLog2FC = c(UpATACLog2FC, DownATACLog2FC))
plotData <- plotData[complete.cases(plotData) & 
                        is.finite(plotData$ExpLog2FC) & 
                        is.finite(plotData$ATACLog2FC), ]

plotData$TypeP <- ifelse(plotData$ExpLog2FC > 0 & plotData$ATACLog2FC > 0, "Up",
                         ifelse(plotData$ExpLog2FC < 0 & plotData$ATACLog2FC < 0, "Down", "Neutral"))

plotData <- plotData[complete.cases(plotData) & 
                        is.finite(plotData$ExpLog2FC) & 
                        is.finite(plotData$ATACLog2FC), ]

fit <- lm(ATACLog2FC ~ ExpLog2FC, data = plotData)

p_value <- summary(fit)$coefficients[2, 4]
r_squared <- summary(fit)$r.squared

p_text <- paste0("P = ", formatC(p_value, format = "e", digits = 2))
r2_text <- paste0("R² = ", round(r_squared, 3))

pdf(args[12], width=5, height=4)

sp <- ggplot(data=plotData, aes(x=ExpLog2FC, y=ATACLog2FC)) +
   geom_smooth(method = "lm", color="steelblue") +
   geom_point(aes(color=TypeP)) +
   scale_color_manual(values=c("Up" = "#FF3333", "Down" = "#A5CCE5", "Neutral" = "grey")) +
   geom_hline(yintercept=0, linetype="dashed", color="grey") +
   geom_vline(xintercept=0, linetype="dashed", color="grey") +
   annotate("text", x=0, y=2, label=r2_text, hjust=0, size=4) +
   annotate("text", x=0, y=1.8, label=p_text, hjust=0, size=4) +
   xlim(-3, 3) +
   ylim(-2, 2) +
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = "black"))

print(sp)
dev.off()