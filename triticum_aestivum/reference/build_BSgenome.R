args <- commandArgs(trailingOnly =TRUE)

twobit_file <- normalizePath(args[1], mustWork=TRUE)
outdir <- args[2]

species_abbr <- "Taestivum"  # e.g., Taestivum
species_name <- ifelse(length(args)>=4, args[4], "Triticum aestivum")
common_name <- ifelse(length(args)>=5, args[5], "wheat")
provider <- ifelse(length(args)>=6, args[6], "IWGSC")
version <- ifelse(length(args)>=7, args[7], "v1")

# Prepare seqs directory
seqdir <- file.path(outdir, "seqs")
dir.create(seqdir, showWarnings=FALSE, recursive=TRUE)
file.copy(twobit_file, file.path(seqdir, basename(twobit_file)), overwrite=TRUE)

library(BSgenome)
library(BSgenomeForge)

# Create seed file
seed_file <- file.path(outdir, paste0("BSgenome.", species_abbr, ".", provider, ".", version, "-seed"))
pkg_name <- paste0("BSgenome.", species_abbr, ".", provider, ".", version)

seed_content <- c(
   paste("Package:", pkg_name),
   paste("Title: Full genome sequences for", species_name),
   paste("Description: Full genome sequences for", species_name, "as provided by", provider, "and stored in Biostrings objects."),
   "Version: 1.0.0",
   paste("organism:", species_name),
   paste("common_name:", common_name),
   paste("provider:", provider),
   paste("provider_version:", version),
   "release_date: 2025",
   "source_url: ftp://your.source.url/",
   paste("organism_biocview:", gsub(" ", "_", species_name)),
   paste("BSgenomeObjname:", species_abbr),
   "seqnames: NULL",
   "circ_seqs: character(0)",
   paste("seqfile_name:", basename(twobit_file))
)

writeLines(seed_content, con=seed_file)
cat("Seed file created:", seed_file, "\n")
file.copy(twobit_file, file.path(target_dir, "single_sequences.2bit"), overwrite = TRUE)
# Build BSgenome package
if (exists("forgeBSgenomeDataPkg", where="package:BSgenomeForge", mode="function")) {
   forgeBSgenomeDataPkg(x=seed_file, seqs_srcdir=seqdir, destdir=outdir, verbose=TRUE)
} else if (exists("forgeBSgenomeDataPkg", where="package:BSgenome", mode="function")) {
   BSgenome::forgeBSgenomeDataPkg(x=seed_file, seqs_srcdir=seqdir, destdir=outdir, verbose=TRUE)
} else {
   stop("forgeBSgenomeDataPkg function not found. Please install BSgenomeForge.")
}
