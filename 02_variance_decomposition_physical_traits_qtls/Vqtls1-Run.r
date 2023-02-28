
rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/2-全基因组方差剖分-2")

for(i in 1:3){

	### STEP1. Load the normalized phenotypes and values of significant QTLs (beta-A, beta-D, MAF) ###
	phedata = read.table(c("1-F2pigs/1_Pigs_NorPhe.txt","2-HSrats/2_HSrats_NorPhe.txt","3-HSmice/3_HSmice_NorPhe.txt")[i], header=T)
	plotdata = read.table(c("1-F2pigs/1-Pigs-SigQTLs.txt","2-HSrats/2-Rats-SigQTLs.txt","3-HSmice/3-Mice-SigQTLs.txt")[i], header=T)

	### STEP2. Calculate Vp using var() based on the normalized phenotypes and compare them with the Vp from GCTA ###
	Vp_norphe = cbind(as.data.frame(sort(unique(plotdata[,1]))), "Vp_nor"=NA); colnames(Vp_norphe)[1] = "Traits"
	for(j in 1:nrow(Vp_norphe)){
		try({  Vp_norphe[j,2] = var(phedata[!is.na(phedata[,Vp_norphe[j,1]]), Vp_norphe[j,1]])  })
	}

	# Vp_GCTA = read.table(c("1-F2pigs/1-Pigs-VpGCTA.txt","2-HSrats/2-Rats-VpGCTA.txt","3-HSmice/3-Mice-VpGCTA.txt")[i], header=T)
	# Vp_norphe = cbind(Vp_norphe, "Vp_GCTA"=NA)
	# for(j in 1:nrow(Vp_norphe)){
	# 	try({  Vp_norphe[j,3] = subset(Vp_GCTA, Vp_GCTA[,1]==Vp_norphe[j,1])[1,2]  })
	# }	 # GCTA uses raw phenotypes to calculate Vp, which is different from the Vp of normalized phenotypes. As the beta values from ADDO was calculated from the normalized phenotypes, so we should Vp_norphe #

	### STEP3. Calculate the Add-Vqtl and Dom-Vqtl for each QTL ###
	plotdata = cbind(plotdata, "Add_Vqtl"=NA, "Dom_Vqtl"=NA)
	for(j in 1:nrow(plotdata)){
		plotdata[j, "Add_Vqtl"] = 2 * plotdata[j,"MAF"] * (1-plotdata[j,"MAF"]) * plotdata[j,"beta_A"] * plotdata[j,"beta_A"]
		plotdata[j, "Dom_Vqtl"] = (2 * plotdata[j,"MAF"] * (1-plotdata[j,"MAF"]) * plotdata[j,"beta_D"])^2
	}

	### STEP4. Calculate the accumulated Add-Vqtl and Dom-Vqtl for each trait ###
	result_tmp = as.data.frame(sort(unique(plotdata[,1])))
	result_tmp = cbind(result_tmp, "Add_Vqtl"=NA, "Dom_Vqtl"=NA, "Add_Var"=NA, "Add_SE"=NA, "Dom_Var"=NA, "Dom_SE"=NA)
	colnames(result_tmp)[1] = "Traits"

	for(k in 1:nrow(result_tmp)){
		plotdata_tmp = subset(plotdata, plotdata[,1]==result_tmp[k,1])
		result_tmp[k,2] = sum(plotdata_tmp[,"Add_Vqtl"])
		result_tmp[k,3] = sum(plotdata_tmp[,"Dom_Vqtl"])
		result_tmp[k,4:7] = as.numeric(plotdata_tmp[1,c("Add_Var", "Add_SE", "Dom_Var", "Dom_SE")])
	}

	### STEP5. Calculate the proportions of Add-Vqtl and Dom-Vqtl to Vp, in case of Vp is not around 1 ###
	res_tmp = cbind(result_tmp[,1], "Vqtla_Vp"=NA, "Vqtld_Vp"=NA, result_tmp[,-1])
	colnames(res_tmp)[1] = "Traits"
	
	for(l in 1:nrow(res_tmp)){
		#try({
				trait_tmp = res_tmp[l, "Traits"]
				res_tmp[l, "Vqtla_Vp"] = res_tmp[l, "Add_Vqtl"]/Vp_norphe[Vp_norphe[,1]==trait_tmp, "Vp_nor"]
				res_tmp[l, "Vqtld_Vp"] = res_tmp[l, "Dom_Vqtl"]/Vp_norphe[Vp_norphe[,1]==trait_tmp, "Vp_nor"]
		#	})
	}

	write.table(res_tmp, file=c("1-F2pigs/1-Pigs-Vqtl.txt","2-HSrats/2-Rats-Vqtl.txt","3-HSmice/3-Mice-Vqtl.txt")[i], row.names=F, col.names=T, quote=F)

}


