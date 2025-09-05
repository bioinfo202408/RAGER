library(tidyr)
library(dplyr)
library(genekitr)
library(Matrix)
library(igraph)
library(dplyr)
library(scales)
args <- commandArgs(trailingOnly = TRUE)

enrich_result <- read.csv(file = args[1], header = TRUE)
enrich_result <- enrich_result[,c(4,7,13)]
enrich_result <- enrich_result %>%
   separate_rows(core_enrichment, sep = "/")
enrich_result <- as.data.frame(enrich_result)
entrez_ids <- enrich_result[, 3]
geneIDs <- transId(
   id = entrez_ids,
   transTo = "ensembl", org = "human", keepNA = FALSE
)
geneIDs <- geneIDs[,2]
expData <- read.csv(args[2], row.names = 1)
exprIDs <- rownames(expData)
exprIDs <- unlist(lapply(strsplit(exprIDs,"\\|"),function(x) x[[1]]))
exprIDs <- gsub("\\.\\d+","",exprIDs)
rownames(expData) <- exprIDs
matchIndexes <- match(geneIDs,rownames(expData))
expData <- expData[matchIndexes,]
expData <- expData[, -1]
entrezIDs <- transId(
   id = rownames(expData),
   transTo = "entrez", org = "human", keepNA = FALSE
)
entrezIDs <- entrezIDs[!duplicated(entrezIDs$input_id), ]
expData$input_id <- rownames(expData)
mergedData <- merge(entrezIDs, expData, by = "input_id")
enrich_result_merged <- enrich_result %>%
   left_join(mergedData %>% dplyr::select(entrezid, log2FoldChange), 
             by = c("core_enrichment" = "entrezid"))
symbolIDs <- transId(
   id = entrez_ids,
   transTo = "symbol", org = "human", keepNA = FALSE
)
symbolIDs <- symbolIDs[!duplicated(symbolIDs$input_id), ]
enrich_result_merged <- merge(
   enrich_result_merged,
   symbolIDs,
   by.x = colnames(enrich_result_merged)[3],
   by.y = colnames(symbolIDs)[1],
   all.x = TRUE
)
gene_set.names <- unique(enrich_result_merged$symbol)
pathway.names <- unique(enrich_result_merged$Description)
node.class <- c(rep("GENE_SET", length(gene_set.names)), 
                rep("PATHWAY", length(pathway.names)))
gene_set.fc <- enrich_result_merged %>%
   group_by(symbol) %>%
   summarize(mean_abs_fc = mean(abs(log2FoldChange))) %>%
   ungroup()
pathway.nes <- enrich_result_merged %>%
   group_by(Description) %>%
   summarize(abs_nes = mean(abs(NES))) %>%
   ungroup()
gene_set_color_pal <- colorRampPalette(c("#99CCFF", "#0066CC", "#003399"))
gene_set_colors <- gene_set_color_pal(100)[as.numeric(cut(
   gene_set.fc$mean_abs_fc, 
   breaks = 100, 
   include.lowest = TRUE
))]
pathway_colors <- rep("#F19A73", length(pathway.names))
node.colors <- c(gene_set_colors, pathway_colors)
gene_set_sizes <- rep(15, length(gene_set.names))
pathway_sizes <- rescale(pathway.nes$abs_nes, to = c(20, 40))
node.sizes <- c(gene_set_sizes, pathway_sizes)
row.indexes <- match(enrich_result_merged$symbol, gene_set.names)
col.indexes <- match(enrich_result_merged$Description, pathway.names) + length(gene_set.names)
enrich.sp.mat <- sparseMatrix(
   i = row.indexes,
   j = col.indexes,
   dims = c(length(c(gene_set.names, pathway.names)), 
            length(c(gene_set.names, pathway.names))),
   dimnames = list(c(gene_set.names, pathway.names),
                   c(gene_set.names, pathway.names))
)
enrich.net.graph <- graph.adjacency(
   enrich.sp.mat,
   mode = "directed"
)
enrich.net.graph <- enrich.net.graph %>%
   set_vertex_attr(name = "class", value = node.class) %>%
   set_vertex_attr(name = "color", value = node.colors) %>%
   set_vertex_attr(name = "size", value = node.sizes) %>%
   set_vertex_attr(name = "log2FC", 
                   value = c(gene_set.fc$mean_abs_fc, 
                             rep(NA, length(pathway.names)))) %>%
   set_vertex_attr(name = "NES", 
                   value = c(rep(NA, length(gene_set.names)), 
                             pathway.nes$abs_nes)) %>%
   set_vertex_attr(name = "shape", 
                   value = ifelse(node.class == "GENE_SET", "ellipse", "rectangle")) %>%
   set_vertex_attr(name = "label", 
                   value = c(gene_set.names, pathway.names))
enrich.net.graph <- set_edge_attr(enrich.net.graph, name = "arrow.size", value = 0.5)
write_graph(enrich.net.graph, args[3], format = "graphml")
