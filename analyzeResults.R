args <- commandArgs(trailingOnly = TRUE)
output_path<-args[1]
plate_meta_path<-args[2]

#output_path<-"/mnt/d/PostdocUnsync/03_crispr/202309_Nelli/crispr202309_data"
#plate_meta_path<-"/mnt/d/PostdocUnsync/03_crispr/202309_Nelli/crispr202309_data/metadata/plateAnnotCrispresso.tsv"

setwd(output_path)

library(ggplot2)
library(reshape2)

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

plateData<-fastRead(plate_meta_path)
plates<-rownames(plateData)

stats<-res<-list();for(P in plates){
    statsTmp<-fastRead(paste0("CRISPAGGR/CRISPRessoAggregate_on_",P,"/CRISPRessoAggregate_mapping_statistics.txt"),row.names = NULL)
    resTmp<-fastRead(paste0("CRISPAGGR/CRISPRessoAggregate_on_",P,"/CRISPRessoAggregate_quantification_of_editing_frequency_by_amplicon.txt"),row.names = NULL)
    sampleNames<-substr(basename(statsTmp$Folder),15,nchar(basename(statsTmp$Folder)))
    statsTmp$sample<-substr(basename(statsTmp$Folder),15,nchar(basename(statsTmp$Folder)))
    resTmp$sample<-substr(basename(resTmp$Folder),15,nchar(basename(resTmp$Folder)))
    
    resTmp$plate<-P
    statsTmp$plate<-P
    
    stats[[P]]<-statsTmp
    res[[P]]<-resTmp
    
};stats<-do.call("rbind",stats);res<-do.call("rbind",res)

resFormated <- data.frame(
    row.names = unique(res$sample),
    plate = res[res$Amplicon == "Reference", "plate"],
    nTotal = res[res$Amplicon == "Reference", "Reads_in_input"],
    nAligned = res[res$Amplicon == "Reference", "Reads_aligned_all_amplicons"]
)
resFormated$nUnaligned <- resFormated$nTotal - resFormated$nAligned

amplicons <- unique(res$Amplicon)
for(ampl in amplicons){
    resOfAmp<-res[res$Amplicon==ampl,]
    rownames(resOfAmp) <- resOfAmp$sample
    nReadCol <- paste0("n",ampl)
    nUnmodifCol <- paste0("unmodified",ampl)
    nModifCol <- paste0("modified",ampl)
    
    resFormated[,nReadCol] <- NA
    resFormated[rownames(resOfAmp),nReadCol]<-resOfAmp$Reads_aligned
    
    resFormated[,nUnmodifCol] <- NA
    resFormated[rownames(resOfAmp),nUnmodifCol]<-resOfAmp$Unmodified
    
    resFormated[,nModifCol] <- NA
    resFormated[rownames(resOfAmp),nModifCol]<-resOfAmp$Modified
}

resFormated$nAmbiguous<- resFormated$nAligned - rowSums(resFormated[,paste0("n",amplicons)],na.rm = TRUE)

fastWrite(res,"results/resultsAggr.tsv",row.names = F)
fastWrite(stats,"results/statsAggr.tsv",row.names = F)
fastWrite(resFormated,"results/resultsFormatted.tsv")

resFormated2plot<-resFormated[resFormated$nAligned>10,c("nUnaligned","nAmbiguous",grep("modif",colnames(resFormated),value=TRUE),"plate")]
resFormated2plot$sample<-rownames(resFormated2plot)

ggData<-reshape2::melt(resFormated2plot,value.name = "n",variable.name="amplicon")

if (nlevels(ggData$amplicon) == 6){
    colorScale <- c("black","grey50","#66B32E","#B8CCA8","#E52421","#E58887")
}else if (nlevels(ggData$amplicon) == 8){
    colorScale <- c("black","grey50","#66B32E","#B8CCA8","#0057A5","#ABD7FF","#E52421","#E58887")
}else{
    colorScale <- rainbow(nlevels(ggData$amplicon))
}


pdf("results/propAmpliconPerPlate.pdf",width = 20,height = 15)
print(ggplot(ggData,mapping=aes(x=sample,y=n,fill=amplicon))+
          geom_bar(position="fill", stat="identity",color="black")+theme_bw()+
          theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold"))+
          facet_wrap(vars(plate),scales = "free_x",ncol = 2)+
          scale_fill_manual(values = colorScale)+
          ylab("%")
)
print(ggplot(ggData,mapping=aes(x=sample,y=n,fill=amplicon))+
          geom_bar(stat="identity",color="black")+theme_bw()+
          theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold"))+
          facet_wrap(vars(plate),scales = "free_x",ncol = 2)+
          scale_fill_manual(values = colorScale)
)
dev.off()
