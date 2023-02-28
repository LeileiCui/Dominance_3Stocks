
rm(list=ls());options(stringsAsFactors=FALSE)
setwd("C:\\Users\\Administrator\\Desktop\\3-HSmice")

ploteQTL = function(dat_AD,dat_A,plotname){

	library(venn)
	
	dat_AD = cbind(dat_AD, "QTL_Name" = 1:nrow(dat_AD))
	dat_AD_1 = subset(dat_AD, dat_AD[,"logP_AD"]>dat_AD[,"thrsuggest"])[,"QTL_Name"]
	dat_AD_2 = subset(dat_AD, dat_AD[,"logP_A"]>dat_AD[,"thrsuggest"])[,"QTL_Name"]
	dat_AD_g = subset(dat_AD, dat_AD[,"logP_AD"]>dat_AD[,"thrgenome"]) # Get the dataset with logP_AD > thrgenome #
	dat_AD_1_g = subset(dat_AD_g, dat_AD_g[,"logP_AD"]>dat_AD_g[,"thrgenome"])[,"QTL_Name"]
	dat_AD_2_g = subset(dat_AD_g, dat_AD_g[,"logP_A"]>dat_AD_g[,"thrgenome"])[,"QTL_Name"]	
	
	dat_A = cbind(dat_A, "QTL_Name" = 1:nrow(dat_A))
	dat_A_1 = subset(dat_A, dat_A[,"logP_AD"]>dat_A[,"thrsuggest"])[,"QTL_Name"]
	dat_A_2 = subset(dat_A, dat_A[,"logP_A"]>dat_A[,"thrsuggest"])[,"QTL_Name"]
	dat_A_g = subset(dat_A, dat_A[,"logP_A"]>dat_A[,"thrgenome"]) # Get the dataset with logP_A > thrgenome #
	dat_A_1_g = subset(dat_A_g, dat_A_g[,"logP_AD"]>dat_A_g[,"thrgenome"])[,"QTL_Name"]
	dat_A_2_g = subset(dat_A_g, dat_A_g[,"logP_A"]>dat_A_g[,"thrgenome"])[,"QTL_Name"]
	
	png(paste0(plotname,"_Venn1.png"),height=6000,width=6000,res=300,type="cairo")
	layout(matrix(1:4,nrow=2),heights=rep(1,2),widths=rep(1,2))
	
	venn_data = list(dat_AD_1,dat_AD_2)
	venn(venn_data, snames=c("AD_Model","A_Model"), ilab=TRUE, zcolor="style", cexil=3, cexsn=3)
	par(mar=c(10,10,10,10))
	title(paste0(plotname," QTLs by AD Model (>",round(mean(dat_AD[,"thrsuggest"]),2),")"),cex.main=3)

	venn_data = list(dat_AD_1_g,dat_AD_2_g)
	venn(venn_data, snames=c("AD_Model","A_Model"), ilab=TRUE, zcolor="style", cexil=3, cexsn=3)
	par(mar=c(10,10,10,10))
	title(paste0(plotname," QTLs by AD Model (>",round(mean(dat_AD[,"thrgenome"]),2),")"),cex.main=3)
	
	venn_data = list(dat_A_1,dat_A_2)
	venn(venn_data, snames=c("AD_Model","A_Model"), ilab=TRUE, zcolor="style", cexil=3, cexsn=3)
	par(mar=c(10,10,10,10))
	title(paste0(plotname," QTLs by A Model (>",round(mean(dat_A[,"thrsuggest"]),2),")"),cex.main=3)

	venn_data = list(dat_A_1_g,dat_A_2_g)
	venn(venn_data, snames=c("AD_Model","A_Model"), ilab=TRUE, zcolor="style", cexil=3, cexsn=3)
	par(mar=c(10,10,10,10))
	title(paste0(plotname," QTLs by A Model (>",round(mean(dat_A[,"thrgenome"]),2),")"),cex.main=3)	
	
	dev.off()

}

QTL_data_AD = read.table("HSmice_QTL_AD.txt",header=T)
QTL_data_A = read.table("HSmice_QTL_A.txt",header=T)
ploteQTL(dat_AD=QTL_data_AD, dat_A=QTL_data_A, plotname="HSmice")

###############################################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)
setwd("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/3-3-显著QTL检测优势对比-1/3-HSmice")

ploteQTL = function(dat_AD,dat_A,plotname){

	library(venn)

	dat_AD = cbind(dat_AD, "QTL_Name" = 1:nrow(dat_AD))
	dat_AD_1 = subset(dat_AD, dat_AD[,"logP_AD"]>dat_AD[,"thrsuggest"])[,"QTL_Name"]
	dat_AD_2 = subset(dat_AD, dat_AD[,"logP_A"]>dat_AD[,"thrsuggest"])[,"QTL_Name"]
	
	dat_A = cbind(dat_A, "QTL_Name" = 1:nrow(dat_A))
	dat_A_1 = subset(dat_A, dat_A[,"logP_AD"]>dat_A[,"thrsuggest"])[,"QTL_Name"]
	dat_A_2 = subset(dat_A, dat_A[,"logP_A"]>dat_A[,"thrsuggest"])[,"QTL_Name"]
	
	png(paste0(plotname,"_Venn2.png"),height=6000,width=6000,res=300,type="cairo")
	layout(matrix(1:4,nrow=2),heights=rep(1,2),widths=rep(1,2))
	
	venn_data = list(dat_AD_1,dat_AD_2)
	venn(venn_data, snames=c("AD_Model","A_Model"), ilab=TRUE, zcolor="style", cexil=3, cexsn=3)
	par(mar=c(10,10,10,10))
	title(paste0(plotname," QTLs by AD Model (>",round(mean(dat_A[,"thrsuggest"]),2),")"),cex.main=3)	
	
	par(mai=c(1.5,1.5,0.5,0.4),mgp=c(5,2,0))
	h <- hist(dat_AD[,"logP_AD"]-dat_AD[,"logP_A"],main="",col=c("grey65"),breaks=50,xlab="logP(AD) - logP(A)",ylab="Number",cex.axis=2.8,cex.lab=3.5)
	box()
	
	venn_data = list(dat_A_1,dat_A_2)
	venn(venn_data, snames=c("AD_Model","A_Model"), ilab=TRUE, zcolor="style", cexil=3, cexsn=3)
	par(mar=c(10,10,10,10))
	title(paste0(plotname," QTLs by A Model (>",round(mean(dat_A[,"thrsuggest"]),2),")"),cex.main=3)

	par(mai=c(1.5,1.5,0.5,0.4),mgp=c(5,2,0))
	h <- hist(dat_A[,"logP_AD"]-dat_A[,"logP_A"],main="",col=c("grey65"),breaks=50,xlab="logP(AD) - logP(A)",ylab="Number",cex.axis=2.8,cex.lab=3.5)
	box()
	
	dev.off()

}

QTL_data_AD = read.table("HSmice_QTL_AD.txt",header=T)
QTL_data_A = read.table("HSmice_QTL_A.txt",header=T)
ploteQTL(dat_AD=QTL_data_AD, dat_A=QTL_data_A, plotname="HSmice1")