library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

args <- commandArgs(trailingOnly = TRUE)
bamfile <- args[1]

bamfile.labels <- gsub(".bam", "", basename(bamfile))
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

outPath <- args[2]
dir.create(outPath)



seqlev <- "chr1" 
seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
which <- as(seqinformation[seqlev], "GRanges")
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
shiftedBamfile <- file.path(outPath, "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
pt <- PTscore(gal1, txs)
pdf(args[3])

plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()

nfr <- NFRscore(gal1, txs)
pdf(args[4])

plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()

tsse <- TSSEscore(gal1, txs)
tsse$TSSEscore
pdf(args[5])

plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()

library(phastCons100way.UCSC.hg38)
txs <- txs[seqnames(txs) %in% "chr1"]
genome <- Hsapiens



objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath,
                              conservation=phastCons100way.UCSC.hg38)

dir(outPath)
objs <- splitGAlignmentsByCut(gal1, txs=txs, outPath = outPath)

library(ChIPpeakAnno)
bamfiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))


pdf(args[6])

cumulativePercentage(bamfiles[1:2], as(seqinformation["chr1"], "GRanges"))
dev.off()

TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)

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

sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))

pdf(args[7])

featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)
dev.off()

pdf(args[8])

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
