library(enrichplot)
args <- commandArgs(trailingOnly = TRUE)
res <- readRDS(args[1])
gsea_result <- res@result
target_pathway_id <- args[2]
pathway_title <- gsea_result$Description[gsea_result$ID == target_pathway_id]
pdf(args[3], width = 10, height = 8)
gseaplot2(
   res,
   geneSetID = target_pathway_id,
   title = pathway_title
)
dev.off()
#Rscript --vanilla GSEA_plot.R gseGO_results.rds GO:0042127 GSEA_plot_for_GO_0042127.pdf