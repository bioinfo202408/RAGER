library(ggplot2)
library(dplyr)
library(ggrepel)
library(cols4all)
args <- commandArgs(trailingOnly = TRUE)
expData <- read.table(file = args[1],
                      sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

exprIDs <- rownames(expData)
exprIDs <- unlist(lapply(strsplit(exprIDs, "\\|"), function(x) x[[1]]))
exprIDs <- gsub("\\.\\d+", "", exprIDs)
rownames(expData) <- exprIDs
expData <- na.omit(expData)
expData <- expData[order(expData$padj), ]

interested_genes_df <- read.table(args[2],
                                  header = FALSE, stringsAsFactors = FALSE)

second_genes_df <- read.table(args[3],
                              header = FALSE, stringsAsFactors = FALSE)

interested_genes_df <- rbind(interested_genes_df, second_genes_df)
library(AnnotationDbi)
library(org.Hs.eg.db)

geneID <- interested_genes_df$V1

g_symbol <- mapIds(org.Hs.eg.db,
                   keys = geneID,
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = "first")

g_symbol <- as.character(g_symbol)
table(is.na(g_symbol))
interested_genes_df$gene_symbol <- g_symbol

colnames(interested_genes_df) <- c("geneID", "geneName")

expData <- expData %>%
   mutate(Significance = case_when(
      rownames(expData) %in% interested_genes_df$geneID ~ "shared_genes",
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Not Significant"
   ))


color_mapping <- c("Up" = "red", "Down" = "blue", "Not Significant" = "grey", "shared_genes" = "green")

expData_filtered <- expData %>%
   filter(abs(log2FoldChange) < 10 & -log10(padj) < 50)

interested_genes_to_label_up <- expData_filtered %>%
   filter(rownames(expData_filtered) %in% interested_genes_df$geneID) %>%
   filter(log2FoldChange > 0) %>%
   arrange(padj) %>%
   head(10)

interested_genes_to_label_down <- expData_filtered %>%
   filter(rownames(expData_filtered) %in% interested_genes_df$geneID) %>%
   filter(log2FoldChange < 0) %>%
   arrange(padj) %>%
   head(10)

interested_genes_to_label <- bind_rows(interested_genes_to_label_up, interested_genes_to_label_down)

interested_genes_to_label$geneID <- rownames(interested_genes_to_label)
interested_genes_to_label <- left_join(interested_genes_to_label, interested_genes_df, by = "geneID")

log2FC_cutoff <- as.numeric(args[4])
padj_cutoff <- as.numeric(args[5])

shared_genes_df <- expData_filtered %>% 
   filter(Significance == "shared_genes")

pdf(args[6], width = 10, height = 7.5)

ggplot() +
   geom_point(data = expData_filtered,
              aes(x = log2FoldChange, y = -log10(padj)),
              color = "grey80", alpha = 0.6, size = 2) +
   
   geom_point(data = shared_genes_df,
              aes(x = log2FoldChange,
                  y = -log10(padj),
                  color = log2FoldChange,
                  size = -log10(pvalue)),
              alpha = 0.8) +
   
   geom_text_repel(data = interested_genes_to_label,
                   aes(x = log2FoldChange, y = -log10(padj), label = geneName),
                   color = "black", size = 6, fontface = "bold", max.overlaps = 100) +
   
   scale_color_continuous_c4a_seq(palette = 'YlGnBu', reverse = TRUE) +
   scale_size_continuous(range = c(1, 6)) +labs(color = "Log2 Fold Change", size = "-Log10(padj)") +
   
   scale_x_continuous(limits = c(-5, 5), breaks = seq(-10, 10, by = 5)) +
   scale_y_continuous(expand = expansion(add = c(2, 0)), limits = c(0, 55), breaks = seq(0, 80, by = 20)) +
   
   geom_hline(yintercept = -log10(padj_cutoff), color = "grey30", linetype = "dotdash", size = 0.6) +
   geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), color = "grey30", linetype = "dashed", size = 0.6) +
   
   labs(title = "Volcano Plot (Shared Genes Highlighted)",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value") +
   theme_minimal() +
   theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 20, hjust = 0.5),
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black")
   )

dev.off()