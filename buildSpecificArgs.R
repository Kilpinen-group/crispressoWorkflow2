args <- commandArgs(trailingOnly = TRUE)
output_path<-args[1]
plate_meta_path<-args[2]
setwd(output_path)

fastRead <- function(fileName, sep = '\t',row.names = 1,as.matrix=FALSE,stringsAsFactors=FALSE,...){
	require(data.table)
	dat <- as.data.frame(data.table::fread(fileName,stringsAsFactors=stringsAsFactors, sep = sep,...))
	if(!is.null(row.names)){
	  rownames(dat) <- dat[,row.names]
	  dat <- dat[,-row.names,drop=FALSE]
	}
	if(as.matrix) dat<-as.matrix(dat)
	return(dat)
}

fastWrite <- function(x, fileName = "default.tsv", headRow="Name",row.names=TRUE,col.names=TRUE, dec=".",sep="\t",...) {
	require(data.table)
	if(is.null(rownames(x))) row.names<-FALSE
	if(is.null(colnames(x))) col.names<-FALSE
	
	if(row.names){
		x=cbind(rownames(x),x)
		colnames(x)[1]<-headRow
	}
	fwrite(x=data.frame(x),file=fileName,sep=sep, row.names = FALSE, col.names = col.names, quote = FALSE, dec=dec,...)
}


dir.create("specificArgs",showWarnings = F)

plateData<-fastRead(plate_meta_path)
sampleNames<-scan("metadata/sampleNames.txt",what="character")

platesPerSample<-substr(sampleNames,1,1)
dataPerSample<-plateData[platesPerSample,]

dataPerSample$sampleName<-sampleNames
rownames(dataPerSample)<-sampleNames

for(sample in sampleNames){
	sampleDat = dataPerSample[sample,]
	cutadaptArg<- paste0("-g ^",sampleDat$adaptFW, " -G ^",sampleDat$adaptRV)
	crispressoArg <- paste0("-a ",sampleDat$original, " -e ",sampleDat$afterHDR, " -fg ",sampleDat$gRNA)
	write(cutadaptArg,paste0("specificArgs/",sample,"_cutadapt.txt"))
	write(crispressoArg,paste0("specificArgs/",sample,"_crispresso.txt"))
}