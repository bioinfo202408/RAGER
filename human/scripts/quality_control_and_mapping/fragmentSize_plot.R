library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
args <- commandArgs(trailingOnly = TRUE)

bamfile <- args[1]
bamfile.labels <- gsub(".bam", "", basename(bamfile))

pdf(args[2])
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()