
rm(list=ls());options(stringsAsFactors=FALSE)

load("Venn_thrsuggest_data1.RData"); load("Venn_thrsuggest_data2.RData")
chr_rm = c(23); phe_rm = c("f41_avgnum","k88ab_avgnum","k88ac_avgnum","k88ad_avgnum","HCT_120")

### Revalue the "MVN_logP_Minor" of SNP with tAB_AA*tAB_BB<0 by NA ###
SNP_list1 = rownames(result_all_all[(!is.na(result_all_all[,"tAB_AA"])) && (!is.na(result_all_all[,"tAB_BB"])),]) # Get the SNP_list1 for SNP with NoNA in "tAB_AA" and "tAB_BB" #
SNP_list2 = rownames(subset(result_all_all[SNP_list1,],result_all_all[SNP_list1,"tAB_AA"]*result_all_all[SNP_list1,"tAB_BB"]<0)) # Based on SNP_list1, generate SNP_list2 for SNP with tAB_AA*tAB_BB<0 #
result_all_all[SNP_list2,"MVN_logP_Minor"] = NA 

### Remove the Sig-SNPs from X chromosome ###
result_A_all = subset(result_A_all,!(result_A_all[,"trait"] %in% phe_rm) & !(result_A_all[,"chr"] %in% chr_rm))
result_D_all = subset(result_D_all,!(result_D_all[,"trait"] %in% phe_rm) & !(result_D_all[,"chr"] %in% chr_rm))
result_T_all = subset(result_T_all,!(result_T_all[,"trait"] %in% phe_rm) & !(result_T_all[,"chr"] %in% chr_rm))
#result_T_all = subset(result_T_all,!(result_T_all[,"trait"] %in% phe_rm) & !(result_T_all[,"chr"] %in% chr_rm) & result_T_all[,"MVN_logP_Minor"]>4.5)
res_A_all = rownames(result_A_all); res_D_all =rownames(result_D_all); res_T_all = rownames(result_T_all)
res_all = sort(unique(c(res_A_all,res_D_all,res_T_all))); names(res_all) = 1:length(res_all)
res_A_allcode = names(res_all[res_all %in% res_A_all])
res_D_allcode = names(res_all[res_all %in% res_D_all])
res_T_allcode = names(res_all[res_all %in% res_T_all])

### Plot1: Using P-values from DvsAD and AvsAD ###
try({

	library(venn); library(vioplot)

	png("ADT-VennVivo_Plot.png",height=3000,width=6000,res=300,type="cairo")
	layout(matrix(c(1,2),ncol=2),widths=c(0.9,1),heights=c(1,1))

	par(mar=c(4,2,5,2), cex.axis=2.5, tck=0.01, mgp=c(3,2,0))
	venn_data = list(res_A_allcode,res_D_allcode,res_T_allcode)
	venn(venn_data, snames=c("A","D","T"), ilab=TRUE, zcolor="style", cexil=2.4, cexsn=2.7, bty="n")
	par(mar=c(4,2,6,2)); title("Overlaps of Sig-SNPs from A/D/T Models",cex.main=2)
	
	result1 = result_A_all; result2 = result_D_all;  result3 = result_T_all
	Pname1 = "X.logP.AddDomCode_Dom."; Pname2 = "X.logP.AddDomCode_Add."
	col1 = "blue"; col2 = "red"
	label1 = "Unique in D"; label2 = "Unique in T"

	SNPset_overlap12 = intersect(rownames(result1),rownames(result2))
	SNPset_overlap13 = intersect(rownames(result1),rownames(result3))
	SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap12]
	SNPset3 = rownames(result3)[!rownames(result3) %in% SNPset_overlap13]

	logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
	logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]
	logP1_set3 = result_all_all[SNPset3,Pname1]; logP1_set3 = logP1_set3[!is.na(logP1_set3)]
	logP2_set3 = result_all_all[SNPset3,Pname2]; logP2_set3 = logP2_set3[!is.na(logP2_set3)]
	
	par(mar=c(4.5,2.5,5,2), cex.axis=2.5, tck=0.01, mgp=c(3,1,0))
	legend1 = c('-log(P-value) of "DvsAD Model"','-log(P-value) of "AvsAD Model"')
	y_lim = max(logP1_set2, logP2_set2, logP1_set3, logP2_set3)*1.2
	vioplot(logP1_set2, logP1_set3, names=c(" "," "), at=c(1,4), col=col1, ylim=c(0,y_lim))
	vioplot(logP2_set2, logP2_set3, names=c(" "," "), at=c(2,5), col=col2, add=TRUE, ylim=c(0,y_lim))
	title("Comparison of -log(P): D and T to A",cex.main=2.2)
	legend("top",inset=.05,legend1,lwd=13, col=c(col1,col2), cex=1.5, horiz=TRUE)
	axis(side=1,at=c(2,5),labels=c(" "," ")); par(mgp=c(3,1.5,0), tck=0); axis(side=1, las=0, at=c(1.5,4.5),labels=c(label1,label2),cex.axis=2)
	dev.off()

})

### Plot2: Using P-values from A and D ###
try({

	library(venn); library(vioplot)

	png("ADT-VennVivo_Plot1.png",height=3000,width=6000,res=300,type="cairo")
	layout(matrix(c(1,2),ncol=2),widths=c(0.9,1),heights=c(1,1))

	par(mar=c(4,2,5,2), cex.axis=2.5, tck=0.01, mgp=c(3,2,0))
	venn_data = list(res_A_allcode,res_D_allcode,res_T_allcode)
	venn(venn_data, snames=c("A","D","T"), ilab=TRUE, zcolor="style", cexil=2.4, cexsn=2.7, bty="n")
	par(mar=c(4,2,6,2)); title("Overlaps of Sig-SNPs from A/D/T Models",cex.main=2)
	
	result1 = result_A_all; result2 = result_D_all;  result3 = result_T_all
	Pname1 = "X.logP.AddCode."; Pname2 = "X.logP.DomCode."
	col1 = "blue"; col2 = "red"
	label1 = "Unique in D"; label2 = "Unique in T"

	SNPset_overlap12 = intersect(rownames(result1),rownames(result2))
	SNPset_overlap13 = intersect(rownames(result1),rownames(result3))
	SNPset2 = rownames(result2)[!rownames(result2) %in% SNPset_overlap12]
	SNPset3 = rownames(result3)[!rownames(result3) %in% SNPset_overlap13]

	logP1_set2 = result_all_all[SNPset2,Pname1]; logP1_set2 = logP1_set2[!is.na(logP1_set2)]
	logP2_set2 = result_all_all[SNPset2,Pname2]; logP2_set2 = logP2_set2[!is.na(logP2_set2)]
	logP1_set3 = result_all_all[SNPset3,Pname1]; logP1_set3 = logP1_set3[!is.na(logP1_set3)]
	logP2_set3 = result_all_all[SNPset3,Pname2]; logP2_set3 = logP2_set3[!is.na(logP2_set3)]
	
	par(mar=c(4.5,2.5,5,2), cex.axis=2.5, tck=0.01, mgp=c(3,1,0))
	legend1 = c('-log(P-value) of "A Model"','-log(P-value) of "D Model"')
	y_lim = max(logP1_set2, logP2_set2, logP1_set3, logP2_set3)*1.2
	vioplot(logP1_set2, logP1_set3, names=c(" "," "), at=c(1,4), col=col1, ylim=c(0,y_lim))
	vioplot(logP2_set2, logP2_set3, names=c(" "," "), at=c(2,5), col=col2, add=TRUE, ylim=c(0,y_lim))
	title("Comparison of -log(P): D and T to A",cex.main=2.2)
	legend("top",inset=.05,legend1,lwd=13, col=c(col1,col2), cex=1.5, horiz=TRUE)
	axis(side=1,at=c(2,5),labels=c(" "," ")); par(mgp=c(3,1.5,0), tck=0); axis(side=1, las=0, at=c(1.5,4.5),labels=c(label1,label2),cex.axis=2)
	dev.off()

})

