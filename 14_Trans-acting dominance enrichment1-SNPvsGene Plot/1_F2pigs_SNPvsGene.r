
# Plot1: Scatter plot of eQTL results

rm(list=ls()); options(stringsAsFactors=FALSE)
setwd("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/6-1-Dom Enrichment-SNPvsGene Plot/1-F2pigs")

plotcols = c("blue3", "deepskyblue", "magenta", "red3") # the color from plotRegion_ddinfo.r #
# plotcols = c("#3C5488FF", "#4DBBD5FF", "orchid3", "#E64B35FF")

### [1] Define the plot function ###
ploteQTL = function(dat,color,plotname){
  
  SNP_chr = as.character(dat$chr);SNP_chr[SNP_chr=="23"] = "19"
  SNP_chr = as.numeric(SNP_chr); SNP_chr[SNP_chr == 0] = NA 
  DGE_chr = dat$DGE_chr; DGE_chr[DGE_chr == "X"] = "19"
  DGE_chr[!(DGE_chr %in% as.character(1:19))] = NA 
  DGE_chr = as.numeric(DGE_chr)
  
  SNP_genomePos = chrmkr[SNP_chr]+dat$pos
  DGE_genomePos = chrmkr[DGE_chr]+dat$DGE_start

  if(color == "red"){ col = rgb(1,0,0,alpha=0.5) }
  if(color == "green"){ col = rgb(0,1,0,alpha=0.5) } 
  if(color == "blue"){ col = rgb(0,0,1,alpha=0.5) }

  if(color == "black"){ col = rgb(66/255,66/255,66/255,alpha=0.5) }
  if(color == "blue3"){ col = rgb(0,0,205/255,alpha=0.5) }
  if(color == "deepskyblue"){ col = rgb(0,191/255,255/255,alpha=0.5) } 
  if(color == "magenta"){ col = rgb(255/255,0,255/255,alpha=0.5) }
  if(color == "red3"){ col = rgb(205/255,0,0,alpha=0.5) }

  # scatter plot eQTL
  png(paste(plotname,"_eQTL",".png",sep=""),height=8850, width=8600,res=600)
  par(mai=c(1.3,1.4,1.5,0.5),mgp=c(3.5,1.5,0))
  plot(SNP_genomePos, DGE_genomePos, pch=19, col = col,
       xaxt="n",yaxt="n",cex.lab=2.5,cex.main=2.9,
	   xlab="SNP genomic positions (CHR)",ylab="Genomic positions of transcripts (CHR)",
	   main=paste0(plotname," :  Genome–genome plot of peak eQTLs"),
       xaxs="i",yaxs="i",cex.axis=2.5,cex=sqrt(dat$logP_AD/5))
  abline(h=chrmkr,lty=2)
  abline(v=chrmkr,lty=2)
  axis(side=1,at=(chrmkr[-1]-(chrLen/2)),labels=c(1:19),cex.axis=2)
  axis(side=2,at=(chrmkr[-1]-(chrLen/2)),labels=c(1:19),cex.axis=2)
  dev.off()

}

### [2] Based on the bim file, generate the length and the cumulative length information for each chr ###
chrpos = read.table("F2pigs.bim",header=F)
idx2keep = !(chrpos[,1] %in% c(0,24))
chr = chrpos[idx2keep,1]; pos = chrpos[idx2keep,4]
chr[chr == 23] = 19; chr_num = 18

chrLen = tapply(pos,chr,max,na.rm=T)
chrmkr = c(0,sapply(1:length(chrLen),function(x){sum(as.numeric(chrLen[1:x]))}))

### [3] Start to plot for each tissue ###
eQTL_data1 = read.table("F2pigs_liver_eQTL_AD.txt",header=T)
eQTL_data1 = eQTL_data1[complete.cases(eQTL_data1),]
eQTL_data1 = subset(eQTL_data1, eQTL_data1[,"DGE_chr"] %in% c(1:chr_num) & eQTL_data1[,"chr"] %in% c(1:chr_num) )

ploteQTL(eQTL_data1,color="black",plotname="F2pigs_liver_AD")
ploteQTL(subset(eQTL_data1,eQTL_data1[,"QTL_Type"]=="additive"),color=plotcols[1],plotname="F2pigs_liver_A")
ploteQTL(subset(eQTL_data1,eQTL_data1[,"QTL_Type"]=="partial-dominant"),color=plotcols[2],plotname="F2pigs_liver_PD")
ploteQTL(subset(eQTL_data1,eQTL_data1[,"QTL_Type"]=="dominant"),color=plotcols[3],plotname="F2pigs_liver_D")
ploteQTL(subset(eQTL_data1,eQTL_data1[,"QTL_Type"]=="overdominant"),color=plotcols[4],plotname="F2pigs_liver_OD")

eQTL_data2 = read.table("F2pigs_muscle_eQTL_AD.txt",header=T)
eQTL_data2 = eQTL_data2[complete.cases(eQTL_data2),]
eQTL_data2 = subset(eQTL_data2, eQTL_data2[,"DGE_chr"] %in% c(1:chr_num) & eQTL_data2[,"chr"] %in% c(1:chr_num) )

ploteQTL(eQTL_data2,color="red",plotname="F2pigs_muscle_AD")
ploteQTL(subset(eQTL_data2,eQTL_data2[,"QTL_Type"]=="additive"),color=plotcols[1],plotname="F2pigs_muscle_A")
ploteQTL(subset(eQTL_data2,eQTL_data2[,"QTL_Type"]=="partial-dominant"),color=plotcols[2],plotname="F2pigs_muscle_PD")
ploteQTL(subset(eQTL_data2,eQTL_data2[,"QTL_Type"]=="dominant"),color=plotcols[3],plotname="F2pigs_muscle_D")
ploteQTL(subset(eQTL_data2,eQTL_data2[,"QTL_Type"]=="overdominant"),color=plotcols[4],plotname="F2pigs_muscle_OD")



