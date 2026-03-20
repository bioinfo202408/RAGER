args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
   cat("Usage: Rscript build_TxDb_OrgDb.R <gtf_file> <outdir>\n")
   quit(status = 1)
}
gtf_file <- normalizePath(args[1], mustWork = TRUE)
outdir <- args[2]
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

library(GenomicFeatures)
library(rtracklayer)
library(AnnotationForge)

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
txdb_file <- file.path(outdir, "wheat_txdb.sqlite")
saveDb(txdb, txdb_file)

gtf <- import(gtf_file)
genes <- gtf[gtf$type == "gene"]
gene_id <- mcols(genes)$gene_id
gene_name <- mcols(genes)$gene_name
if (is.null(gene_name)) gene_name <- gene_id
gene_info <- data.frame(
   GID = as.character(gene_id),
   SYMBOL = as.character(gene_name),
   GENENAME = as.character(gene_name),
   stringsAsFactors = FALSE
)
gene_info <- gene_info[!is.na(gene_info$GID) & gene_info$GID != "", ]
gene_info$SYMBOL[is.na(gene_info$SYMBOL) | gene_info$SYMBOL == ""] <- gene_info$GID[is.na(gene_info$SYMBOL) | gene_info$SYMBOL == ""]
gene_info$GENENAME[is.na(gene_info$GENENAME) | gene_info$GENENAME == ""] <- gene_info$GID[is.na(gene_info$GENENAME) | gene_info$GENENAME == ""]
gene_info <- gene_info[!duplicated(gene_info$GID), ]

go_df <- data.frame(
   GID = character(),
   GO = character(),
   EVIDENCE = character(),
   stringsAsFactors = FALSE
)

makeOrgPackage(
   gene_info = gene_info,
   go = go_df,
   version = "1.0",
   maintainer = "user <user@email.com>",
   author = "user",
   outputDir = outdir,
   tax_id = "4565",
   genus = "Triticum",
   species = "aestivum",
   goTable = "go"
)

cat("TxDb and OrgDb built successfully.\n")