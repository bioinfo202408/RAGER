library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(txdbmaker)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("Usage: Rscript script.R <gtf_file> <output_bed> <dnase_rep1> <dnase_rep2> <k27ac_rep1> <k27ac_rep2> <k4me1_rep1> <k4me1_rep2> <k4me3>")
}

gtf_file <- args[1]
output_file <- args[2]
dnase1_file <- args[3]
dnase2_file <- args[4]
k27ac1_file <- args[5]
k27ac2_file <- args[6]
k4me11_file <- args[7]
k4me12_file <- args[8]
k4me3_file <- args[9]

read_bed <- function(bed_path) {
  bed <- import.bed(bed_path)
  return(bed)
}

dnase1 <- read_bed(dnase1_file)
dnase2 <- read_bed(dnase2_file)
k27ac1 <- read_bed(k27ac1_file)
k27ac2 <- read_bed(k27ac2_file)
k4me11 <- read_bed(k4me11_file)
k4me12 <- read_bed(k4me12_file)
k4me3 <- read_bed(k4me3_file)

dnase_consensus   <- union(dnase1, dnase2)
cat("Consensus DNase peaks:", length(dnase_consensus), "\n")

k27ac_union <- union(k27ac1, k27ac2)
cat("Union H3K27ac peaks:", length(k27ac_union), "\n")

k4me1_union <- union(k4me11, k4me12)
cat("Union H3K4me1 peaks:", length(k4me1_union), "\n")

markers_union <- union(k27ac_union, k4me1_union)
candidate_enhancer <- intersect(dnase_consensus, markers_union)
cat("Candidate enhancer peaks:", length(candidate_enhancer), "\n")

txdb <- makeTxDbFromGFF(
   file = gtf_file,
   format = "gtf",
   organism = "Triticum aestivum"
)

tss <- promoters(txdb, upstream = 0, downstream = 1)
promoter_region <- promoters(txdb, upstream = 1000, downstream = 0)
promoter_region <- reduce(promoter_region)
cat("Promoter regions:", length(promoter_region), "\n")

enhancer_no_promoter <- setdiff(candidate_enhancer, promoter_region)
final_enhancer <- setdiff(enhancer_no_promoter, k4me3)
cat("Final enhancer-like regions:", length(final_enhancer), "\n")

final_enhancer$name <- paste0("wheat_enhancer_", seq_along(final_enhancer))
final_enhancer$score <- 100

export.bed(final_enhancer, output_file)
cat("Analysis completed. Output file:", output_file, "\n")

#Rscript make_enhancer.R \
#    /path/to/Triticum_aestivum.IWGSC.62.gtf \
#    /path/to/output/wheat_enhancer_like_annotation.bed \
#    /path/to/GSM3564342_macs_CS_seedlings_DNaseI_seq_rep1_peaks.bed \
#    /path/to/GSM3564343_macs_CS_seedlings_DNaseI_seq_rep2_peaks.bed \
#    /path/to/GSM3449726_macs_CS_seedlings_H3K27ac_ChIP_seq_rep1_peaks.bed \
#    /path/to/GSM3449727_macs_CS_seedlings_H3K27ac_ChIP_seq_rep2_peaks.bed \
#    /path/to/GSM3449724_macs_CS_seedlings_H3K4me1_ChIP_seq_rep1_peaks.bed \
#    /path/to/GSM3449725_macs_CS_seedlings_H3K4me1_ChIP_seq_rep2_peaks.bed \
#    /path/to/GSM3449722_macs_CS_seedlings_H3K4me3_ChIP_seq_peaks.bed