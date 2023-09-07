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
	statsTmp$sample<-sampleNames
	resTmp$sample<-sampleNames
	
	resTmp$plate<-P
	statsTmp$plate<-P
	
	stats[[P]]<-statsTmp
	res[[P]]<-resTmp
	
};stats<-do.call("rbind",stats);res<-do.call("rbind",res)

resFormated<-data.frame(row.names = unique(res$sample))
resFormated$nAligned<-res[res$Amplicon=="Reference","Reads_aligned_all_amplicons"]
resFormated$plate<-res[res$Amplicon=="Reference","plate"]

for(col in unique(res$Amplicon)){
	resFormated[,paste0("n",col)]<-res[res$Amplicon==col,"Reads_aligned"]
	resFormated[,paste0("unmodified",col)]<-res[res$Amplicon==col,"Unmodified"]
	resFormated[,paste0("modified",col)]<-res[res$Amplicon==col,"Modified"]
}

fastWrite(res,"results/resultsAggr.tsv",row.names = F)
fastWrite(stats,"results/statsAggr.tsv",row.names = F)
fastWrite(resFormated,"results/resultsFormatted.tsv")

resFormated2plot<-resFormated[resFormated$nAligned>10,c(grep("modif",colnames(resFormated),value=TRUE),"plate")]
resFormated2plot$sample<-rownames(resFormated2plot)

ggData<-reshape2::melt(resFormated2plot,value.name = "n",variable.name="seqType")

pdf("results/propAmpliconPerPlate.pdf",width = 20,height = 15)
print(ggplot(ggData,mapping=aes(x=sample,y=n,fill=seqType))+
	geom_bar(position="fill", stat="identity",color="black")+theme_bw()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold"))+
	facet_wrap(vars(plate),scales = "free_x",ncol = 2)+
	scale_fill_manual(values = c("#E52421","#E58887","#66B32E","#B8CCA8"))
)
dev.off()
