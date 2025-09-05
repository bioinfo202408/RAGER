library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genekitr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
Shared_UPenhancer <- read.table(args[1],header=FALSE,stringsAsFactors=FALSE)
enh_anno <- read.table(args[2],header=FALSE,stringsAsFactors=FALSE)
colnames(Shared_UPenhancer) <- c("EnhancerID")
colnames(enh_anno) <- c("chr", "start", "end","EnhancerID")

Shared_UPenhancer_nano <- merge(Shared_UPenhancer, enh_anno[, c("EnhancerID", "chr", "start", "end")], 
                                by = "EnhancerID", all.x = TRUE)
Shared_UPenhancer_nano <- na.omit(Shared_UPenhancer_nano)
gr <- GRanges(
   seqnames = Shared_UPenhancer_nano$chr,
   ranges = IRanges(Shared_UPenhancer_nano$start, Shared_UPenhancer_nano$end)
)

extracted_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
names(extracted_seqs) <- Shared_UPenhancer_nano$EnhancerID
writeXStringSet(extracted_seqs, args[3])

Shared_DOWNenhancer <- read.table(args[4],header=FALSE,stringsAsFactors=FALSE)

colnames(Shared_DOWNenhancer) <- c("EnhancerID")
colnames(enh_anno) <- c("chr", "start", "end","EnhancerID")

Shared_DOWNenhancer_nano <- merge(Shared_DOWNenhancer, enh_anno[, c("EnhancerID", "chr", "start", "end")], 
                                by = "EnhancerID", all.x = TRUE)
Shared_DOWNenhancer_nano <- na.omit(Shared_DOWNenhancer_nano)
gr <- GRanges(
   seqnames = Shared_DOWNenhancer_nano$chr,
   ranges = IRanges(Shared_DOWNenhancer_nano$start, Shared_DOWNenhancer_nano$end)
)

extracted_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
names(extracted_seqs) <- Shared_DOWNenhancer_nano$EnhancerID
writeXStringSet(extracted_seqs, args[5])
