
rm(list=ls());options(stringsAsFactors=FALSE)
setwd("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/6-2-Dom Enrichment-Venn Plot/3-HSmice")

ploteQTL_2Venn = function(dat,plotname){

	library(venn)
	
	dat = cbind(dat, "QTL_Name" = 1:nrow(dat))
	Add_eQTLs = subset(dat, dat[,1] %in% c("additive","partial-dominant"))[,"QTL_Name"]
	Dom_eQTLs = subset(dat, dat[,1] %in% c("dominant","overdominant"))[,"QTL_Name"]
	Cis_eQTLs =  subset(dat, dat[,3] == "cis")[,"QTL_Name"]
	Trans_eQTLs = subset(dat, dat[,3] == "trans")[,"QTL_Name"]
	
	png(paste0(plotname,"_Venn.png"),height=3000,width=6000,res=300,type="cairo")
	layout(matrix(1:2,nrow=1),heights=2,widths=rep(1,2))
	
	venn_data = list(Add_eQTLs,Cis_eQTLs)
	venn(venn_data, snames=c("Additive-eQTLs","Cis-eQTLs"), ilab=F, zcolor="style", ilcs=2.5, sncs=2)

	par(mar=c(10,10,10,10))
	title(paste0(plotname," eQTLs"),cex.main=2.2)
	
	venn_data = list(Dom_eQTLs,Trans_eQTLs)
	venn(venn_data, snames=c("Dominant-eQTLs","Trans-eQTLs"), ilab=F, zcolor="style", ilcs=2.5, sncs=2)

	par(mar=c(10,10,10,10))
	title(paste0(plotname," eQTLs"),cex.main=2.2)

	dev.off()

}

ploteQTL_4Venn = function(dat,plotname){

	library(venn); library(vioplot)
	
	dat = cbind(dat, "QTL_Name" = 1:nrow(dat))
	Add_eQTLs = subset(dat, dat[,1] %in% c("additive","partial-dominant"))[,"QTL_Name"]
	Dom_eQTLs = subset(dat, dat[,1] %in% c("dominant","overdominant"))[,"QTL_Name"]
	Cis_eQTLs =  subset(dat, dat[,3] == "cis")[,"QTL_Name"]
	Trans_eQTLs = subset(dat, dat[,3] == "trans")[,"QTL_Name"]
	
	venn_data = list(Add_eQTLs,Dom_eQTLs,Cis_eQTLs,Trans_eQTLs)
	
	png(paste0(plotname,"_Venn.png"),height=3000,width=3000,res=300)
	venn(venn_data, snames=c("Add","Dom","Cis","Trans"), ilab=F, zcolor="style", ilcs=2.7, sncs=3.2)
	
	par(mar=c(10,10,10,10))
	title(paste0(plotname," eQTLs"),cex.main=3.3)
	
	dev.off()

}

eQTL_data1 = read.table("HSmice_hippocampus_eQTL_AD.txt",header=T)
eQTL_data1 = eQTL_data1[complete.cases(eQTL_data1),]
eQTL_data1 = cbind(eQTL_data1, "AddDom_Code"=0, "CisTrans_Code"=0)
eQTL_data1[eQTL_data1[,"eQTL_Type"]%in%c("partial-dominant", "overdominant"), "AddDom_Code"]=1
eQTL_data1[eQTL_data1[,"cis_trans"]=="trans", "CisTrans_Code"]=1
eQTL_data1_plot = eQTL_data1[,c("eQTL_Type", "AddDom_Code", "cis_trans", "CisTrans_Code")]

ploteQTL_2Venn(eQTL_data1_plot,plotname="HSmice_hippocampus_AvsD")
ploteQTL_4Venn(eQTL_data1_plot,plotname="HSmice_hippocampus")

eQTL_data2 = read.table("HSmice_liver_eQTL_AD.txt",header=T)
eQTL_data2 = eQTL_data2[complete.cases(eQTL_data2),]
eQTL_data2 = cbind(eQTL_data2, "AddDom_Code"=0, "CisTrans_Code"=0)
eQTL_data2[eQTL_data2[,"eQTL_Type"]%in%c("partial-dominant", "overdominant"), "AddDom_Code"]=1
eQTL_data2[eQTL_data2[,"cis_trans"]=="trans", "CisTrans_Code"]=1
eQTL_data2_plot = eQTL_data2[,c("eQTL_Type", "AddDom_Code", "cis_trans", "CisTrans_Code")]

ploteQTL_2Venn(eQTL_data2_plot,plotname="HSmice_liver_AvsD")
ploteQTL_4Venn(eQTL_data2_plot,plotname="HSmice_liver")

eQTL_data3 = read.table("HSmice_lung_eQTL_AD.txt",header=T)
eQTL_data3 = eQTL_data3[complete.cases(eQTL_data3),]
eQTL_data3 = cbind(eQTL_data3, "AddDom_Code"=0, "CisTrans_Code"=0)
eQTL_data3[eQTL_data3[,"eQTL_Type"]%in%c("partial-dominant", "overdominant"), "AddDom_Code"]=1
eQTL_data3[eQTL_data3[,"cis_trans"]=="trans", "CisTrans_Code"]=1
eQTL_data3_plot = eQTL_data3[,c("eQTL_Type", "AddDom_Code", "cis_trans", "CisTrans_Code")]

ploteQTL_2Venn(eQTL_data3_plot,plotname="HSmice_lung_AvsD")
ploteQTL_4Venn(eQTL_data3_plot,plotname="HSmice_lung")

fisher.test(eQTL_data1_plot[,2], eQTL_data1_plot[,4])$p.value
chisq.test(eQTL_data1_plot[,2], eQTL_data1_plot[,4])$p.value

fisher.test(eQTL_data2_plot[,2], eQTL_data2_plot[,4])$p.value
chisq.test(eQTL_data2_plot[,2], eQTL_data2_plot[,4])$p.value

fisher.test(eQTL_data3_plot[,2], eQTL_data3_plot[,4])$p.value
chisq.test(eQTL_data3_plot[,2], eQTL_data3_plot[,4])$p.value


