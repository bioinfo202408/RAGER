library(Matrix)
library(igraph)
library(viridis)
library(org.Mm.eg.db)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
target_tfs <- unlist(strsplit(args[3], ","))

tf_motifEnrich <- read.table(file=input_file, 
                             sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="")
tf_motifEnrich <- subset(tf_motifEnrich, V5 == "Yes")

tf_motifEnrich_filtered <- tf_motifEnrich %>%
   filter(V2 %in% target_tfs)

filtered_top10_network <- tf_motifEnrich_filtered %>%
   group_by(V2) %>%
   arrange(desc(V4)) %>%
   slice_head(n = 10) %>%
   ungroup()

geneID <- filtered_top10_network$V3 %>% 
   as.character() %>% 
   na.omit()

g_symbol <- mget(geneID, org.Mm.egSYMBOL, ifnotfound=NA) %>% 
   as.character()

filtered_top10_network$V3 <- g_symbol
tf_motifEnrich_filtered <- filtered_top10_network[, c(2, 3, 4, 6)] %>% 
   as.data.frame()
colnames(tf_motifEnrich_filtered) <- c("TF", "Gene", "Enrichment_Score", "P_value")

tf.names <- unique(tf_motifEnrich_filtered$TF)
gene.names <- unique(tf_motifEnrich_filtered$Gene)
node.class <- c(rep("TF", length(tf.names)), 
                rep("GENE", length(gene.names)))

tf.pvalues <- aggregate(P_value ~ TF, 
                        data = tf_motifEnrich_filtered, 
                        FUN = min)

generate_color <- function(pvals) {
   neg_log_p <- -log10(pvals)
   color_palette <- colorRampPalette(c("#FF9999", "#FF6666", "#990000"))
   norm_p <- (neg_log_p - min(neg_log_p)) / (max(neg_log_p) - min(neg_log_p))
   color_palette(100)[cut(norm_p, breaks = 100, include.lowest = TRUE)]
}

tf_colors <- generate_color(tf.pvalues$P_value)
node.colors <- rep("#AAD7F4", length(c(tf.names, gene.names))) 
node.colors[1:length(tf.names)] <- tf_colors[match(tf.names, tf.pvalues$TF)]

gene.sizes <- aggregate(Enrichment_Score ~ Gene, 
                        data = tf_motifEnrich_filtered, 
                        FUN = max)

normalize_size <- function(x) {
   min_size <- 5
   max_size <- 20
   (x - min(x)) / (max(x) - min(x)) * (max_size - min_size) + min_size
}

node.sizes <- rep(10, length(c(tf.names, gene.names)))
gene_size_values <- normalize_size(gene.sizes$Enrichment_Score)
gene_positions <- match(gene.sizes$Gene, c(tf.names, gene.names))
node.sizes[gene_positions] <- gene_size_values

row.indexes <- match(tf_motifEnrich_filtered$TF, c(tf.names, gene.names))
col.indexes <- match(tf_motifEnrich_filtered$Gene, c(tf.names, gene.names))

tf.gene.sp.mat <- sparseMatrix(
   i = col.indexes,
   j = row.indexes,
   x = 1,
   symmetric = FALSE,
   dims = c(length(c(tf.names, gene.names)), 
            length(c(tf.names, gene.names))),
   dimnames = list(c(tf.names, gene.names), 
                   c(tf.names, gene.names))
)

tf.gene.net.graph <- graph.adjacency(
   tf.gene.sp.mat,
   mode = "directed",
   weighted = NULL
)

tf.gene.net.graph <- tf.gene.net.graph %>%
   set_vertex_attr(name = "class", value = node.class) %>%
   set_vertex_attr(name = "color", value = node.colors) %>%
   set_vertex_attr(name = "pvalue", 
                   value = c(tf.pvalues$P_value[match(tf.names, tf.pvalues$TF)], 
                             rep(NA, length(gene.names)))) %>%
   set_vertex_attr(name = "shape", 
                   value = ifelse(node.class == "TF", "square", "circle")) %>%
   set_vertex_attr(name = "size", value = node.sizes)

tf.gene.net.graph <- tf.gene.net.graph %>%
   set_edge_attr(name = "arrow.size", value = 0.5)

write_graph(tf.gene.net.graph, 
            output_file, 
            format = "graphml")
