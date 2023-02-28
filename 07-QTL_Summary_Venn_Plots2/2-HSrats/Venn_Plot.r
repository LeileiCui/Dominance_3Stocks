
rm(list=ls());options(stringsAsFactors=FALSE)

Stock_name = c("HSrats"); data_dir = c("/SAN/mottlab/heterosis/3_HSrats"); Sig_level = "suggest"

### Part1: Calculate the sum of significant SNPs in 7 models A/D/AD/AD_A/AD_D/Indi/T and their overlaps SNP num ###
AD_dir = paste0(data_dir,"/3NonAdd_QTL/1_A+D_0.999999r2/2_Pvalue"); T_dir = paste0(data_dir,"/3NonAdd_QTL/2_T-statistics_0.999999r2/2_Pvalue")

AD_files = dir(AD_dir,pattern=".txt.tar.gz"); AD_names = gsub("Pvalue_","",gsub(".txt.tar.gz","",AD_files)) 
T_files = dir(T_dir,pattern=".txt.tar.gz"); T_names = gsub("Pvalue_","",gsub(".txt.tar.gz","",T_files)) 
Map_Info = read.table(paste0(data_dir,"/3NonAdd_QTL/1_A+D_0.999999r2/1_PheGen/0_Data/Data_clean.map"),header=F)

result_A_all = c(); result_D_all = c(); result_AD_all = c(); result_AD_A_all = c(); result_AD_D_all = c(); result_Indi_all = c(); result_T_all = c(); result_all_all = c()
for(i in 1:length(AD_files)){
try({
	
	setwd(AD_dir)
	system(paste0("tar zxvf ",AD_files[i]))
	tmp_AD_raw = read.table(gsub(".tar.gz","",AD_files[i]),header=T); rownames(tmp_AD_raw) = paste(AD_names[i],Map_Info[,2],sep=":")
	tmp_AD = subset(tmp_AD_raw,tmp_AD_raw[,"chr"]!=0 & tmp_AD_raw[,"pos"]!=0) # Remove SNPs without CHR and POS info #
	tmp_AD = tmp_AD[complete.cases(tmp_AD),] # Remove SNPs with NA results #
	system(paste0("rm ",gsub(".tar.gz","",AD_files[i])))

	setwd(T_dir)
	system(paste0("tar zxvf ",T_files[i]))
	tmp_T_raw = read.table(gsub(".tar.gz","",T_files[i]),header=T); rownames(tmp_T_raw) =  paste(AD_names[i],Map_Info[,2],sep=":")
	tmp_T = subset(tmp_T_raw,tmp_T_raw[,"chr"]!=0 & tmp_T_raw[,"pos"]!=0) # Remove SNPs without CHR and POS info #
	tmp_T = tmp_T[complete.cases(tmp_T),] # Remove SNPs with NA results #
	tmp_T = subset(tmp_T, tmp_T[,"tAB_AA"]*tmp_T[,"tAB_BB"]>0) # Remove SNPs with tAB_AA*tAB_BB<0 # 
	system(paste0("rm ",gsub(".tar.gz","",T_files[i])))

	thrgenome_AD = -log10(0.05/nrow(tmp_AD)); thrsuggest_AD = -log10(1/nrow(tmp_AD))
	thrgenome_T = -log10(0.05/nrow(tmp_T)); thrsuggest_T = -log10(1/nrow(tmp_T))
	if(Sig_level == "genome"){
		thrused_AD = thrgenome_AD; thrused_T = thrgenome_T
	}else{
		thrused_AD = thrsuggest_AD; thrused_T = thrsuggest_T
	}
	
	result_A = subset(tmp_AD,tmp_AD[,"X.logP.AddCode."] > thrused_AD)
	result_D = subset(tmp_AD,tmp_AD[,"X.logP.DomCode."] > thrused_AD)
	result_AD = subset(tmp_AD,tmp_AD[,"X.logP.AddDomCode_Null."] > thrused_AD)
	result_AD_A = subset(tmp_AD,tmp_AD[,"X.logP.AddDomCode_Add."] > thrused_AD)
	#result_AD_D = subset(tmp_AD,tmp_AD[,"X.logP.AddDomCode_Dom."] > thrused_AD)
	result_Indi = subset(tmp_AD,tmp_AD[,"X.logP.IndiCode."] > thrused_AD)
	result_T = subset(tmp_T,tmp_T[,"MVN_logP_Minor"] > thrused_T)
	
	if(nrow(result_A)>0){ result_A = cbind("trait" = AD_names[i], result_A); result_A_all = rbind(result_A_all,result_A) }
	if(nrow(result_D)>0){ result_D = cbind("trait" = AD_names[i], result_D); result_D_all = rbind(result_D_all,result_D) }
	if(nrow(result_AD)>0){ result_AD = cbind("trait" = AD_names[i], result_AD); result_AD_all = rbind(result_AD_all,result_AD) }
	if(nrow(result_AD_A)>0){ result_AD_A = cbind("trait" = AD_names[i], result_AD_A); result_AD_A_all = rbind(result_AD_A_all,result_AD_A) }
	#if(nrow(result_AD_D)>0){ result_AD_D = cbind("trait" = AD_names[i], result_AD_D); result_AD_D_all = rbind(result_AD_D_all,result_AD_D) }
	if(nrow(result_Indi)>0){ result_Indi = cbind("trait" = AD_names[i], result_Indi); result_Indi_all = rbind(result_Indi_all,result_Indi) }
	if(nrow(result_T)>0){ result_T = cbind("trait" = T_names[i], result_T); result_T_all = rbind(result_T_all,result_T) }
	
	tmp_AD_T = cbind(tmp_AD_raw,tmp_T_raw[,-(1:2)]); rownames(tmp_AD_T) = paste(AD_names[i],Map_Info[,2],sep=":")
	result_all_all = rbind(result_all_all, tmp_AD_T)
	
	print(i)
})
}

res_A_all = rownames(result_A_all); res_D_all =rownames(result_D_all); res_AD_all =rownames(result_AD_all); res_AD_A_all = rownames(result_AD_A_all); res_AD_D_all = rownames(result_AD_D_all); res_Indi_all = rownames(result_Indi_all); res_T_all = rownames(result_T_all)

res_all = sort(unique(c(res_A_all,res_D_all,res_AD_all,res_AD_A_all,res_AD_D_all,res_Indi_all,res_T_all))); names(res_all) = 1:length(res_all)
res_A_allcode = names(res_all[res_all %in% res_A_all])
res_D_allcode = names(res_all[res_all %in% res_D_all])
res_AD_allcode = names(res_all[res_all %in% res_AD_all])
res_AD_A_allcode = names(res_all[res_all %in% res_AD_A_all])
res_AD_D_allcode = names(res_all[res_all %in% res_AD_D_all])
res_Indi_allcode = names(res_all[res_all %in% res_Indi_all])
res_T_allcode = names(res_all[res_all %in% res_T_all])

Overlap_A = c(AvsD = length(intersect(res_A_all,res_D_all)), AvsAD = length(intersect(res_A_all,res_AD_all)), AvsAD_A = length(intersect(res_A_all,res_AD_A_all)), AvsAD_D = length(intersect(res_A_all,res_AD_D_all)), AvsIndi = length(intersect(res_A_all,res_Indi_all)), AvsT = length(intersect(res_A_all,res_T_all)))
Overlap_D = c(DvsAD = length(intersect(res_D_all,res_AD_all)), DvsAD_A = length(intersect(res_D_all,res_AD_A_all)), DvsAD_D = length(intersect(res_D_all,res_AD_D_all)), DvsIndi = length(intersect(res_D_all,res_Indi_all)), DvsT = length(intersect(res_D_all,res_T_all)))
Overlap_AD = c(ADvsAD_A = length(intersect(res_AD_all,res_AD_A_all)), ADvsAD_D = length(intersect(res_AD_all,res_AD_D_all)), ADvsIndi = length(intersect(res_AD_all,res_Indi_all)), ADvsT = length(intersect(res_AD_all,res_T_all)))
Overlap_AD_A = c(AD_AvsAD_D = length(intersect(res_AD_A_all,res_AD_D_all)), AD_AvsIndi = length(intersect(res_AD_A_all,res_Indi_all)), AD_AvsT = length(intersect(res_AD_A_all,res_T_all)))
Overlap_AD_D = c(AD_DvsIndi = length(intersect(res_AD_D_all,res_Indi_all)), AD_DvsT = length(intersect(res_AD_D_all,res_T_all)))
Overlap_Indi = c(IndivsT = length(intersect(res_Indi_all,res_T_all)))

setwd(paste0(AD_dir,"/../../Venn-Plot"))
save(result_A_all,result_D_all,result_AD_all,result_AD_A_all,result_AD_D_all,result_Indi_all,result_T_all,result_all_all, file=paste0("Venn_thr",Sig_level,"_data1.RData"))
save(res_A_allcode,res_D_allcode,res_AD_allcode,res_AD_A_allcode,res_AD_D_allcode,res_Indi_allcode,res_T_allcode, file=paste0("Venn_thr",Sig_level,"_data2.RData"))

### Part2: Venn Plot ###
library(venn)

png(paste0("All-",Stock_name,"_Venn_thr",Sig_level,".png"),height=3000,width=3000,res=300,type="cairo")
venn_data = list(res_A_allcode,res_AD_allcode,res_D_allcode,res_AD_A_allcode,res_T_allcode)
venn(venn_data, snames=c("A","A+D","D","AvsA+D","T"), ilab=TRUE, zcolor="style", cexil=1.4, cexsn=1.8)
par(mar=c(7,7,7,7))
#title(paste0("Significant SNPs Overlaps in different models (",Stock_name,")"),cex.main=1.7)
dev.off()

### Part2 4in1 Plot #

### Revalue the "MVN_logP_Minor" of SNP with tAB_AA*tAB_BB<0 by NA ###
SNP_list1 = rownames(result_all_all[(!is.na(result_all_all[,"tAB_AA"])) && (!is.na(result_all_all[,"tAB_BB"])),]) # Get the SNP_list1 for SNP with NoNA in "tAB_AA" and "tAB_BB" #
SNP_list2 = rownames(subset(result_all_all[SNP_list1,],result_all_all[SNP_list1,"tAB_AA"]*result_all_all[SNP_list1,"tAB_BB"]<0)) # Based on SNP_list1, generate SNP_list2 for SNP with tAB_AA*tAB_BB<0 #
result_all_all[SNP_list2,"MVN_logP_Minor"] = NA 

### Plot1: A/AD/D ###
try({

	### Part2: Venn Plot ###
	library(venn); library(vioplot)

	png(paste0(Stock_name,"_Venn_thr",Sig_level,"0.png"),height=3000,width=4500,res=300,type="cairo")
	layout(matrix(c(1,1,1,2,3,4),ncol=2),widths=c(1,1),heights=c(1,1,1))

	venn_data = list(res_A_allcode,res_AD_allcode,res_D_allcode)
	venn(venn_data, snames=c("Add","Add+Dom","Dom"), ilab=TRUE, zcolor="style", cexil=2.4, cexsn=2.8)
	#title(paste0("Significant SNPs Overlaps in various models (",Stock_name,")"),cex.main=1.7)

	### Part3-1: Box Plot of -logP value ###
		result1 = result_A_all; result2 = result_AD_all
		Pname1 = "X.logP.AddCode."; Pname2 = "X.logP.AddDomCode_Null."
		col1 = "red1"; col2 = "red4"
		title1 = "-logP of Sig-SNPs Comparison: Add Model with Add+Dom Model"
		legend1 = c("-log(P-value) from Add Model","-log(P-value) from Add+Dom Model")
		label1 = "Sig-SNPs only from Add"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from A+D"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4, bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.5)

	### Part3-2: Box Plot of -logP value ###
		result1 = result_D_all; result2 = result_AD_all
		Pname1 = "X.logP.DomCode."; Pname2 = "X.logP.AddDomCode_Null."
		col1 = "blue1"; col2 = "blue4"
		title1 = "-logP of Sig-SNPs Comparison: Dom Model with Add+Dom Model"
		legend1 = c("-log(P-value) from Dom Model","-log(P-value) from Add+Dom Model")
		label1 = "Sig-SNPs only from Dom"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from A+D"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4, bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.5)
		
	### Part3-3: Box Plot of -logP value ###
		result1 = result_A_all; result2 = result_D_all
		Pname1 = "X.logP.AddCode."; Pname2 = "X.logP.DomCode."
		col1 = "green1"; col2 = "green4"
		title1 = "-logP of Sig-SNPs Comparison: Add Model with Dom Model"
		legend1 = c("-log(P-value) from Add Model","-log(P-value) from Dom Model")
		label1 = "Sig-SNPs only from Add"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from Dom"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4 ,bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.5)
		
	dev.off()

})

### Plot2: T/ADvsA/D ###
try({

	### Part2: Venn Plot ###
	library(venn); library(vioplot)

	png(paste0(Stock_name,"_Venn_thr",Sig_level,"1.png"),height=3000,width=4500,res=300,type="cairo")
	layout(matrix(c(1,1,1,2,3,4),ncol=2),widths=c(1,1),heights=c(1,1,1))

	venn_data = list(res_T_allcode,res_AD_A_allcode,res_D_allcode)
	venn(venn_data, snames=c("T","Add+Dom vs Add","Dom"), ilab=TRUE, zcolor="style", cexil=2.4, cexsn=2.8)
	#title(paste0("Significant SNPs Overlaps in various models (",Stock_name,")"),cex.main=1.7)

	### Part3-1: Box Plot of -logP value ###
		result1 = result_T_all; result2 = result_AD_A_all
		Pname1 = "MVN_logP_Minor"; Pname2 = "X.logP.AddDomCode_Add."
		col1 = "red1"; col2 = "red4"
		title1 = "-logP of Sig-SNPs Comparison: T Model with Add+Dom vs Add Model"
		legend1 = c("-log(P-value) from T Model","-log(P-value) from Add+Dom vs Add Model")
		label1 = "Sig-SNPs only from T"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from ADvsA"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4, bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)

	### Part3-2: Box Plot of -logP value ###
		result1 = result_D_all; result2 = result_AD_A_all
		Pname1 = "X.logP.DomCode."; Pname2 = "X.logP.AddDomCode_Add."
		col1 = "blue1"; col2 = "blue4"
		title1 = "-logP of Sig-SNPs Comparison: Dom Model with Add+Dom vs Add Model"
		legend1 = c("-log(P-value) from Dom Model","-log(P-value) from Add+Dom vs Add Model")
		label1 = "Sig-SNPs only from D"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from ADvsA"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4, bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)
		
	### Part3-3: Box Plot of -logP value ###
		result1 = result_T_all; result2 = result_D_all
		Pname1 = "MVN_logP_Minor"; Pname2 = "X.logP.DomCode."
		col1 = "green1"; col2 = "green4"
		title1 = "-logP of Sig-SNPs Comparison: T Model with Dom Model"
		legend1 = c("-log(P-value) from T Model","-log(P-value) from Dom Model")
		label1 = "Sig-SNPs only from T"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from D"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4 ,bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)
		
	dev.off()

})

### Plot3: ADvsA/AD/T ###
try({

	### Part2: Venn Plot ###
	library(venn); library(vioplot)

	png(paste0(Stock_name,"_Venn_thr",Sig_level,"2.png"),height=3000,width=4500,res=300,type="cairo")
	layout(matrix(c(1,1,1,2,3,4),ncol=2),widths=c(1,1),heights=c(1,1,1))

	venn_data = list(res_AD_allcode,res_AD_A_allcode,res_T_allcode)
	venn(venn_data, snames=c("Add+Dom","Add+Dom vs Add","T"), ilab=TRUE, zcolor="style", cexil=2.4, cexsn=2.7)
	#title(paste0("Significant SNPs Overlaps in various models (",Stock_name,")"),cex.main=1.7)

	### Part3-1: Box Plot of -logP value ###
		result1 = result_AD_all; result2 = result_AD_A_all
		Pname1 = "X.logP.AddCode."; Pname2 = "X.logP.AddDomCode_Add."
		col1 = "red1"; col2 = "red4"
		title1 = "-logP of Sig-SNPs Comparison: Add+Dom Model with Add+Dom vs Add Model"
		legend1 = c("-log(P-value) from Add+Dom Model","-log(P-value) from Add+Dom vs Add Model")
		label1 = "Sig-SNPs only from AD"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from ADvsA"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4, bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)

	### Part3-2: Box Plot of -logP value ###
		result1 = result_AD_all; result2 = result_T_all
		Pname1 = "X.logP.AddCode."; Pname2 = "MVN_logP_Minor"
		col1 = "blue1"; col2 = "blue4"
		title1 = "-logP of Sig-SNPs Comparison: Add+Dom Model with T Model"
		legend1 = c("-log(P-value) from Add+Dom Model","-log(P-value) from T Model")
		label1 = "Sig-SNPs only from AD"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from T"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4, bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)
		
	### Part3-3: Box Plot of -logP value ###
		result1 = result_AD_A_all; result2 = result_T_all
		Pname1 = "X.logP.AddDomCode_Add."; Pname2 = "MVN_logP_Minor"
		col1 = "green1"; col2 = "green4"
		title1 = "-logP of Sig-SNPs Comparison: Add+Dom vs Add Model with T Model"
		legend1 = c("-log(P-value) from Add+Dom vs Add Model","-log(P-value) from T Model")
		label1 = "Sig-SNPs only from ADvsA"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from T"

		SNPset_overlap = intersect(rownames(result1),rownames(result2))
		SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		title(title1,cex.main=1.6)
		legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4 ,bty="n")
		axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)
		
	dev.off()

})

# ### Plot4: ADvsA/AD/ADvsD ###
# try({

	# ### Part2: Venn Plot ###
	# library(venn); library(vioplot)

	# png(paste0(Stock_name,"_Venn_thr",Sig_level,"3.png"),height=3000,width=4500,res=300,type="cairo")
	# layout(matrix(c(1,1,1,2,3,4),ncol=2),widths=c(1,1),heights=c(1,1,1))

	# venn_data = list(res_AD_A_allcode,res_AD_allcode,res_AD_D_allcode)
	# venn(venn_data, snames=c("Add+Dom vs Add","Add+Dom","Add+Dom vs Dom"), ilab=TRUE, zcolor="style", cexil=2.4, cexsn=2.2)
	# #title(paste0("Significant SNPs Overlaps in various models (",Stock_name,")"),cex.main=1.7)

	# ### Part3-1: Box Plot of -logP value ###
		# result1 = result_AD_all; result2 = result_AD_A_all
		# Pname1 = "X.logP.AddCode."; Pname2 = "X.logP.AddDomCode_Add."
		# col1 = "red1"; col2 = "red4"
		# title1 = "-logP of Sig-SNPs Comparison: Add+Dom Model with Add+Dom vs Add Model"
		# legend1 = c("-log(P-value) from Add+Dom Model","-log(P-value) from Add+Dom vs Add Model")
		# label1 = "Sig-SNPs only from AD"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from ADvsA"

		# SNPset_overlap = intersect(rownames(result1),rownames(result2))
		# SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		# SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		# logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		# logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		# logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		# logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		# logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		# logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		# par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		# y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		# vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		# vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		# title(title1,cex.main=1.6)
		# legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4, bty="n")
		# axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)

	# ### Part3-2: Box Plot of -logP value ###
		# result1 = result_AD_all; result2 = result_AD_D_all
		# Pname1 = "X.logP.AddCode."; Pname2 = "X.logP.AddDomCode_Dom."
		# col1 = "blue1"; col2 = "blue4"
		# title1 = "-logP of Sig-SNPs Comparison: Add+Dom Model with Add+Dom vs Dom Model"
		# legend1 = c("-log(P-value) from Add+Dom Model","-log(P-value) from Add+Dom vs Dom Model")
		# label1 = "Sig-SNPs only from AD"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from ADvsD"

		# SNPset_overlap = intersect(rownames(result1),rownames(result2))
		# SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		# SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		# logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		# logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		# logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		# logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		# logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		# logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		# par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		# y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		# vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		# vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		# title(title1,cex.main=1.6)
		# legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.4, bty="n")
		# axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)
		
	# ### Part3-3: Box Plot of -logP value ###
		# result1 = result_AD_A_all; result2 = result_AD_D_all
		# Pname1 = "X.logP.AddDomCode_Add."; Pname2 = "X.logP.AddDomCode_Dom."
		# col1 = "green1"; col2 = "green4"
		# title1 = "-logP of Sig-SNPs Comparison: Add+Dom vs Add Model with  Add+Dom vs Dom Model"
		# legend1 = c("-log(P-value) from Add+Dom vs Add Model","-log(P-value) from Add+Dom vs Dom Model")
		# label1 = "Sig-SNPs only from ADvsA"; label2 = "Sig-SNPs from Both"; label3 = "Sig-SNPs only from ADvsD"

		# SNPset_overlap = intersect(rownames(result1),rownames(result2))
		# SNPset1 = rownames(result1)[!rownames(result1) %in% SNPset_overlap]
		# SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap]

		# logP1_set1 = result_all_all[SNPset1,Pname1]; logP1_set1 = logP1_set1[!is.na(logP1_set1)]
		# logP2_set1 = result_all_all[SNPset1,Pname2]; logP2_set1 = logP2_set1[!is.na(logP2_set1)]
		# logP1_set12 =  result_all_all[SNPset_overlap,Pname1]; logP1_set12 = logP1_set12[!is.na(logP1_set12)]
		# logP2_set12 = result_all_all[SNPset_overlap,Pname2]; logP2_set12 = logP2_set12[!is.na(logP2_set12)]
		# logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
		# logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]

		# par(mar=c(2.5,2.5,2.5,2),cex.lab=2, cex.axis=1.7, tck=0.01)
		# y_lim = max(logP1_set1, logP1_set12, logP1_set2, logP2_set1, logP2_set12, logP2_set2)*1.05
		# vioplot(logP1_set1, logP1_set12, logP1_set2, names=c(" "," "," "), at=c(1,4,7), col=col1, ylim=c(0,y_lim))
		# vioplot(logP2_set1, logP2_set12, logP2_set2, names=c(" "," "," "), at=c(2,5,8), col=col2, add=TRUE, ylim=c(0,y_lim))
		# title(title1,cex.main=1.4)
		# legend("topleft",inset=.05,legend1,lwd=10, col=c(col1,col2), cex=1.2 ,bty="n")
		# axis(side=1, at=c(2,5,8),labels=c(" "," "," ")); par(tck=0); axis(side=1, las=0, at=c(1.5,4.5,7.5),labels=c(label1,label2,label3),cex.axis = 1.4)
		
	# dev.off()

# })

