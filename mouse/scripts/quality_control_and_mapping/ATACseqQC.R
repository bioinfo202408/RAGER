library(ATACseqQC)
library(BSgenome.Mmusculus.UCSC.mm10)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
args <- commandArgs(trailingOnly = TRUE)
bamfile <- args[1]

bamfile.labels <- gsub(".bam", "", basename(bamfile))
pdf(args[2])
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                    "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                  "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                  "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                  "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                  "U2"))
library(Rsamtools)
bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags
#outPath <- "/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263045/splited"
outPath <- args[3]
dir.create(outPath)
## shift the coordinates of 5'ends of alignments in the bam file
## if you don't have an available TxDb, please refer
## GenomicFeatures::makeTxDbFromGFF to create one from gff3 or gtf file.
seqlev <- "chr1" ## subsample data for quick run
seqinformation <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)
which <- as(seqinformation[seqlev], "GRanges")
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
shiftedBamfile <- file.path(outPath, "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
pt <- PTscore(gal1, txs)
pdf(args[4])
#pdf("/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263045/PTscore_plot.pdf")
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()

nfr <- NFRscore(gal1, txs)
pdf(args[5])
#pdf("/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263045/NFRscore_plot.pdf")
plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()

tsse <- TSSEscore(gal1, txs)
tsse$TSSEscore
pdf(args[6])
#pdf("/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263045/TSSEscore_plot.pdf")
plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()

txs <- txs[seqnames(txs) %in% "chr1"]
genome <- Mmusculus

objs <- splitGAlignmentsByCut(gal1, txs=txs, outPath = outPath)

library(ChIPpeakAnno)
bamfiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))
## Plot the cumulative percentage of tag allocation in nucleosome-free 
## and mononucleosome bam files.
pdf(args[7])
#pdf("/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263045/cumulative_percentage_plot.pdf")
cumulativePercentage(bamfiles[1:2], as(seqinformation["chr1"], "GRanges"))
dev.off()

TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
## estimate the library size for normalization
(librarySize <- estLibSize(bamfiles))

NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)
## log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
#plot heatmap
pdf(args[8])
#pdf("/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263045/heatmap.pdf")
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)
dev.off()

pdf(args[9])
#pdf("/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263045/coverage_curve_plot.pdf")
out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()
