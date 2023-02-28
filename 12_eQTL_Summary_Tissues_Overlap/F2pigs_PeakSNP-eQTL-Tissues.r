
##################################################################################################
### Filter overlap QTLs in two methods and SNPs with biggest P-value that are both >thrsuggest ###
##################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

AD_dir = "/SAN/mottlab/heterosis/1_F2pigs/3NonAdd_eQTL/1_A+D_L/SNP_Peak" 
AD_e_dir = "/SAN/mottlab/heterosis/1_F2pigs/3NonAdd_eQTL/1_A+D_M/SNP_Peak"; outname = "liver_muscle" 

# AD_dir = "/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL/1_A+D_0.999999r2_Amygdala/SNP_Peak"
# AD_e_dir = "/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL/1_A+D_0.999999r2_Heart/SNP_Peak"; outname = "amygdala_heart_t"  

# AD_dir = "/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL_Gene/1_A+D_amygdala/SNP_Peak"
# AD_e_dir = "/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL_Gene/1_A+D_heart/SNP_Peak"; outname = "amygdala_heart_g" 

# AD_dir = "/SAN/mottlab/heterosis/4_HSmice/3NonAdd_eQTL_Res/1_A+D_hippo/SNP_Peak"
# AD_e_dir = "/SAN/mottlab/heterosis/4_HSmice/3NonAdd_eQTL_Res/1_A+D_liver/SNP_Peak"; outname = "hippo_liver" 

# AD_dir = "/SAN/mottlab/heterosis/4_HSmice/3NonAdd_eQTL_Res/1_A+D_hippo/SNP_Peak"
# AD_e_dir = "/SAN/mottlab/heterosis/4_HSmice/3NonAdd_eQTL_Res/1_A+D_lung/SNP_Peak"; outname = "hippo_lung" 

# AD_dir = "/SAN/mottlab/heterosis/4_HSmice/3NonAdd_eQTL_Res/1_A+D_liver/SNP_Peak"
# AD_e_dir = "/SAN/mottlab/heterosis/4_HSmice/3NonAdd_eQTL_Res/1_A+D_lung/SNP_Peak"; outname = "liver_lung" 

Distance = 2000000; thres_sig = 6.5
result_all1 = read.table(paste0(AD_dir,"/PeakSNP_AD_thrsuggest.txt"),header=T)
result_all2 = read.table(paste0(AD_e_dir,"/PeakSNP_AD_thrsuggest.txt"),header=T)

result_all1 = subset(result_all1, result_all1[,"logP_AD"]>thres_sig)
result_all2 = subset(result_all2, result_all2[,"logP_AD"]>thres_sig)

xxx = cbind(result_all2[,1:3],do.call(rbind,strsplit(result_all2[,"CI_r2"],"-"))); colnames(xxx) = c("trait","chr","pos","CI_START","CI_END") # Use xxx stands for eQTL results #
xxx[,"CI_START"] = as.numeric(xxx[,"CI_START"]); xxx[,"CI_END"] = as.numeric(xxx[,"CI_END"])
yyy = cbind(result_all1[,1:3],do.call(rbind,strsplit(result_all1[,"CI_r2"],"-"))); colnames(yyy) = c("trait","chr","pos","CI_START","CI_END") # Use yyy stands for QTL results #
yyy[,"CI_START"] = as.numeric(yyy[,"CI_START"]); yyy[,"CI_END"] = as.numeric(yyy[,"CI_END"])

xxx = xxx[!duplicated(paste(xxx[,"trait"],xxx[,"chr"],xxx[,"pos"],sep=":")),]
yyy = yyy[!duplicated(paste(yyy[,"trait"],yyy[,"chr"],yyy[,"pos"],sep=":")),]

# [1] Get the overlap QTLs #
tmp_result_all = data.frame()
for(i in 1:nrow(yyy)){ # Based on QTLs, to search the overlap eQTLs from T-statistics results #

	tmp_xxx = subset(xxx,xxx[,"chr"] == yyy[i,"chr"])
	if(dim(tmp_xxx)[1] == 0){
		tmp_result = NA
	}else{	
		tmp_result = c()
		for(j in 1:nrow(tmp_xxx)){
			if( ((tmp_xxx[j,"CI_START"] >= yyy[i,"CI_START"]-Distance) & (tmp_xxx[j,"CI_START"] <= yyy[i,"CI_START"]+Distance)) | ((tmp_xxx[j,"CI_END"] >= yyy[i,"CI_END"]-Distance) & (tmp_xxx[j,"CI_END"] <= yyy[i,"CI_END"]+Distance)) ){
				tmp_res = cbind(paste(yyy[i,"trait"],yyy[i,"chr"],yyy[i,"pos"],yyy[i,"CI_START"],yyy[i,"CI_END"],sep=":"), paste(tmp_xxx[j,],collapse=":"))
				tmp_result = rbind(tmp_result,tmp_res)
			}
		}
	}
	
	tmp_result_all = rbind(tmp_result_all,tmp_result)
	print(i)

}
colnames(tmp_result_all) = c("PeakSNP_QTL","PeakSNP_eQTL")
result_all = tmp_result_all[complete.cases(tmp_result_all),]

# [2] Get the common significant SNPs (just choose the peak SNP by "A+D Results"/"QTL Results") for each pair #
result_all = cbind(result_all,Common_Peak=NA); result_detail = c()
for(i in 1:nrow(result_all)){

	tmp_qtl1 = strsplit(result_all[i,1],":")[[1]]
	tmp_qtl2 = strsplit(result_all[i,2],":")[[1]]
	names(tmp_qtl1) = names(tmp_qtl2) = c("trait","chr","pos","CI_START","CI_END")
	
	setwd(AD_dir); system(paste0("tar zxvf ../2_Pvalue/Pvalue_",tmp_qtl1[1],".txt.tar.gz"))
	tmp_pvalue1 = read.table(paste0("Pvalue_",tmp_qtl1[1],".txt"),header=T)
	tmp_pvalue1 = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]!=0 & tmp_pvalue1[,"pos"]!=0) # Remove SNPs without CHR and POS info #
	tmp_pvalue1 = tmp_pvalue1[complete.cases(tmp_pvalue1),] # Remove SNPs with NA results #
	thrgenome1 = -log10(0.05/nrow(tmp_pvalue1)); thrsuggest1 = -log10(1/nrow(tmp_pvalue1))
	system(paste0("rm Pvalue_",tmp_qtl1[1],".txt"))
	
	setwd(AD_e_dir); system(paste0("tar zxvf ../2_Pvalue/Pvalue_",tmp_qtl2[1],".txt.tar.gz"))
	tmp_pvalue2 = read.table(paste0("Pvalue_",tmp_qtl2[1],".txt"),header=T)
	tmp_pvalue2 = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]!=0 & tmp_pvalue2[,"pos"]!=0) # Remove SNPs without CHR and POS info #
	tmp_pvalue2 = tmp_pvalue2[complete.cases(tmp_pvalue2),] # Remove SNPs with NA results #
	thrgenome2 = -log10(0.05/nrow(tmp_pvalue2)); thrsuggest2 = -log10(1/nrow(tmp_pvalue2))
	system(paste0("rm Pvalue_",tmp_qtl2[1],".txt"))	
	
	if(tmp_qtl1[3] == tmp_qtl2[3]){ # Case1: If tmp_qtl1 and tmp_qtl2 share the common peak SNP #
		tmp_qtl = tmp_qtl1; result_all[i,3] = result_all[i,1]
		result_detail_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue1[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
		result_detail_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue2[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
		result_detail = rbind(result_detail,cbind(trait=tmp_qtl1[1],result_detail_AD,DGE=tmp_qtl2[1],result_detail_T))
	}else{ # Case2: If tmp_qtl1 and tmp_qtl2 just have overlap CI region #
		region_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl1["chr"]) & tmp_pvalue1[,"pos"]>=as.numeric(tmp_qtl1["CI_START"]) & tmp_pvalue1[,"pos"]<=as.numeric(tmp_qtl1["CI_END"])) # Get the region SNPs for tmp_qtl1 #
		region_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl2["chr"]) & tmp_pvalue2[,"pos"]>=as.numeric(tmp_qtl2["CI_START"]) & tmp_pvalue2[,"pos"]<=as.numeric(tmp_qtl2["CI_END"])) # Get the region SNPs for tmp_qtl2 #
		
		common_SNP = region_AD[,"pos"][region_AD[,"pos"] %in% region_T[,"pos"]]
		region_AD_common = region_AD[region_AD[,"pos"] %in% common_SNP,]
		region_AD_common_Sig = region_AD_common[region_AD_common[,"X.logP.AddDomCode_Null."]>thrsuggest1,] # Get the common region & significant SNPs for tmp_qtl1 #
		region_T_common = region_T[region_T[,"pos"] %in% common_SNP,]
		region_T_common_Sig = region_T_common[region_T_common[,"X.logP.AddDomCode_Null."]>thrsuggest2,] # Get the common region & significant SNPs for tmp_qtl2 #
		
		if(nrow(region_AD_common_Sig) != 0 && nrow(region_T_common_Sig) != 0){
			common_SNP_Sig = region_AD_common_Sig[,"pos"][region_AD_common_Sig[,"pos"] %in% region_T_common_Sig[,"pos"]]
			if(length(common_SNP_Sig) > 0){
				region_AD_common_Key = region_AD_common_Sig[region_AD_common_Sig[,"pos"] %in% common_SNP_Sig,]
				tmp_qtl =c(tmp_qtl1[1], as.character(region_AD_common_Key[region_AD_common_Key[,"X.logP.AddDomCode_Null."]==max(region_AD_common_Key[,"X.logP.AddDomCode_Null."]),1:2][1,]))
				names(tmp_qtl) = c("trait","chr","pos")
				result_all[i,3] = paste(tmp_qtl,collapse=":")
				result_detail_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue1[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
				result_detail_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue2[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
				result_detail = rbind(result_detail,cbind(trait=tmp_qtl1[1],result_detail_AD,DGE=tmp_qtl2[1],result_detail_T))
			}
		}
	}
	print(i)
}	

setwd(paste0(AD_e_dir,"/../../Overlap-Tissues"))
write.table(result_all,file=paste0("PeakSNP-QTL_eQTL-", outname, "-ADCommon1.txt"),row.names=F,col.names=T,quote=F) # Contains all the overlap QTL-eQTL pairs and the 3rd col is the "Both-Sig SNP" #
write.table(result_detail,file=paste0("PeakSNP-QTL_eQTL-", outname, "-ADCommon2.txt"),row.names=F,col.names=T,quote=F) # Just contains the "Both-Sig SNPs" from overlap QTL-eQTL pairs #

