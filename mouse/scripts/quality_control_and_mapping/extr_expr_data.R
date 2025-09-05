args = commandArgs(trailingOnly=TRUE)
pattern <- file.path(args[1], "*", "gene_abund.tab")
files <- Sys.glob(pattern)
expData <- read.delim(file=files[1],header=TRUE,stringsAsFactors=FALSE)
expData <- expData[,c(1,2,9)]
posIndex <- regexpr('stringtie.*gene_abund',files[1])
sampleid <- substring(files[1],posIndex+14,posIndex+attr(posIndex,'match.length')-12)
uniqIndex <- which(!duplicated(as.character(expData[,1])))
uniqEnsemids <- as.character(expData[uniqIndex,1])
expData <- expData[uniqIndex,]
columnNames <- c("gene_id","gene_name",sampleid)
for (filename in files[2:length(files)]){
	if(file.info(filename)[1] > 0){
		sampleMat <- read.delim(file=filename,header=TRUE,stringsAsFactors=FALSE)
		matchIndex <- match(uniqEnsemids,sampleMat[,1])
		expData <- cbind(expData,sampleMat[matchIndex,9])
		posIndex <- regexpr('stringtie.*gene_abund',filename)
		sampleid <- substring(filename,posIndex+14,posIndex+attr(posIndex,'match.length')-12)
		cat(sampleid,"\n")
		columnNames <- append(columnNames,sampleid)
	}
}
colnames(expData) <- columnNames
write.table(expData,file=args[2],sep='\t',quote=FALSE,row.names=T,col.names=T)
