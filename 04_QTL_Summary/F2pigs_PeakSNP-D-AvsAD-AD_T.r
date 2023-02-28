
rm(list=ls());options(stringsAsFactors=FALSE)

chr_rm = c(23); phe_rm = c("f41_avgnum","k88ab_avgnum","k88ac_avgnum","k88ad_avgnum","HCT_120")

AD_dir = "/SAN/mottlab/heterosis/1_F2pigs/3NonAdd_QTL/1_A+D" 
T_dir = "/SAN/mottlab/heterosis/1_F2pigs/3NonAdd_QTL/2_T-statistics" 
folder_name = "1_F2pigs"; file_name = "F2";
CI_Pvalue_decrease = 2; CI_Pvalue_r2 = 0.8; Distance = 2000000 # Distance allow QTL and QTL(or eQTL) are shared one #

##################################################################################################
### Part1: PeakSNP-D_T.r (Based on D P-value to choose the QTLs from A+D Model)
##################################################################################################

try({

	setwd(AD_dir)
	system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --recode --out tmp"))
	system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --freq --out tmp_freq"))
	Map_Info = read.table("tmp.map",header=F); Freq_Info = read.table("tmp_freq.frq",header=T)
	Freq_Info = cbind(Map_Info,Freq_Info)
	system("rm tmp*")

	AD_results = dir(paste0(AD_dir,"/2_Pvalue"),pattern=".txt.tar.gz")
	AD_phes = gsub("Pvalue_","",gsub(".txt.tar.gz","",AD_results)) 
	T_results = dir(paste0(T_dir,"/2_Pvalue"),pattern=".txt.tar.gz")
	T_phes = gsub("Pvalue_","",gsub(".txt.tar.gz","",T_results)) 

	# STEP1: Extract potential QTLs>thrsuggest based on A+D Results #
	setwd(AD_dir); result_all1 = c()
	for(i in 1:length(AD_phes)){
	try({

		system(paste0("tar zxvf 2_Pvalue/",AD_results[i]))
		tmp_AD_raw = read.table(gsub(".tar.gz","",AD_results[i]),header=T)
		tmp_AD = subset(tmp_AD_raw,tmp_AD_raw[,"chr"]!=0 & tmp_AD_raw[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_AD = tmp_AD[complete.cases(tmp_AD),] # Remove SNPs with NA results #
		system(paste0("rm ",gsub(".tar.gz","",AD_results[i])))
		
		system(paste0("tar zxvf ",T_dir,"/2_Pvalue/",AD_results[i]))
		tmp_T_raw = read.table(gsub(".tar.gz","",AD_results[i]),header=T)
		system(paste0("rm ",gsub(".tar.gz","",AD_results[i])))
		
		thrgenome = -log10(0.05/nrow(tmp_AD)); thrsuggest = -log10(1/nrow(tmp_AD))
		
		if(max(tmp_AD[,"X.logP.DomCode."]) > thrsuggest){
			tmp_result_all1 = c()
			tmp_AD_Sig = subset(tmp_AD,tmp_AD[,"X.logP.DomCode."]  > thrsuggest)
			for(j in 1:length(unique(tmp_AD_Sig[,"chr"]))){
				# [1] Get the A+D Model Info for the Peak SNP #
				tmp_AD_Sig_chr = subset(tmp_AD_Sig,tmp_AD_Sig[,"chr"] == unique(tmp_AD_Sig[,"chr"])[j])
				tmp_AD_Sig_max = subset(tmp_AD_Sig_chr,tmp_AD_Sig_chr[,"X.logP.DomCode."] == max(tmp_AD_Sig_chr[,"X.logP.DomCode."]))[1,]
				tmp_AD_result = cbind(AD_phes[i], tmp_AD_Sig_max[,"chr"], tmp_AD_Sig_max[,"pos"], tmp_AD_Sig_max[,c("X.logP.AddDomCode_Add.","X.logP.AddCode.","X.logP.DomCode.","X.logP.AddDomCode_Null.","X.logP.IndiCode.","betaAdd","betaDom")], tmp_AD_Sig_max[,"betaDom"]/tmp_AD_Sig_max[,"betaAdd"])
			
				# [2] Get the T-statistics Model Info for the Peak SNP #
				tmp_T_Sig_max = subset(tmp_T_raw,tmp_T_raw[,"chr"]==tmp_AD_Sig_max[,"chr"] & tmp_T_raw[,"pos"]==tmp_AD_Sig_max[,"pos"])
				tmp_T_result = tmp_T_Sig_max[,c("MVN_logP_Minor","tAB_AA","tAB_BB","betaAB_AA","betaAB_BB","varAB_AA","varAB_BB","df_res")]
			
				# [3] Get the Confidence Level for the Peak SNP #
				tmp_Sig_SNP = subset(Freq_Info,Freq_Info[,1]==tmp_AD_Sig_max[,"chr"] & Freq_Info[,4]==tmp_AD_Sig_max[,"pos"])[,"SNP"]
				system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --r2 --ld-snp ",tmp_Sig_SNP," --ld-window 99999 --ld-window-kb 50000 --ld-window-r2 0"," --out tmp_Sig_SNP"))
				TEMP.ld <- read.table(file="tmp_Sig_SNP.ld",header=T,stringsAsFactors=F); rownames(TEMP.ld) = TEMP.ld[,"SNP_B"]; system("rm tmp_Sig_SNP.*")
				TEMP.ld.tmp = TEMP.ld[TEMP.ld[,"R2"]>CI_Pvalue_r2,]
				tmp_CI_result1 = paste(TEMP.ld.tmp[1,"BP_B"],TEMP.ld.tmp[nrow(TEMP.ld.tmp),"BP_B"],sep="-")
				
				tmp_AD_chr = subset(tmp_AD,tmp_AD[,"chr"] == unique(tmp_AD_Sig[,"chr"])[j] & tmp_AD[,"X.logP.DomCode."] >= (tmp_AD_Sig_max[,"X.logP.DomCode."]-CI_Pvalue_decrease))
				tmp_CI_result2 = paste(tmp_AD_chr[1,"pos"],tmp_AD_chr[nrow(tmp_AD_chr),"pos"],sep="-")
			
				# [4] Get the MAF for the Peak SNP #
				tmp_MAF_result = subset(Freq_Info,Freq_Info[,1]==tmp_AD_Sig_max[,"chr"] & Freq_Info[,4]==tmp_AD_Sig_max[,"pos"])[,"MAF"]
				
				tmp_result = cbind(tmp_AD_result,tmp_T_result,tmp_MAF_result,tmp_CI_result1,tmp_CI_result2,thrgenome,thrsuggest)
				tmp_result_all1 = rbind(tmp_result_all1,tmp_result)
			}
			colnames(tmp_result_all1) = c("trait","chr","pos","logP_AD_A","logP_A","logP_D","logP_AD","logP_Indi","beta_Add","beta_OverDom","OD/A","MVN_logP_Minor","tAB_AA","tAB_BB","betaAB_AA","betaAB_BB","varAB_AA","varAB_BB","df_res","MAF","CI_r2","CI_logP","thrgenome","thrsuggest")
		
			print(paste0("Phe",i," : ",AD_phes[i]," has significant SNP in ",length(unique(tmp_AD_Sig[,"chr"]))," chrs!"))
		}else{
			tmp_result_all1 = c()
			print(paste0("Phe",i," : ",AD_phes[i]," no significant SNP." ))
		}
		
		result_all1 = rbind(result_all1,tmp_result_all1)
		
	})
	}
	result_all1 = subset(result_all1,!(result_all1[,"trait"] %in% phe_rm) & !(result_all1[,"chr"] %in% chr_rm))
	write.table(result_all1,file="PeakSNP_D.txt",row.names=F,col.names=T,quote=F)

	# STEP2: Extract potential QTLs>thrsuggest based on T-statistics Results #
	setwd(T_dir); result_all2 = c()

	if(!file.exists("PeakSNP_T.txt")){

		for(i in 1:length(T_phes)){
		try({

			system(paste0("tar zxvf 2_Pvalue/",T_results[i]))
			tmp_T_raw = read.table(gsub(".tar.gz","",T_results[i]),header=T)
			tmp_T = subset(tmp_T_raw,tmp_T_raw[,"chr"]!=0 & tmp_T_raw[,"pos"]!=0) # Remove SNPs without CHR and POS info #
			tmp_T = tmp_T[complete.cases(tmp_T),] # Remove SNPs with NA results #
			tmp_T = subset(tmp_T, tmp_T[,"tAB_AA"]*tmp_T[,"tAB_BB"]>0)	
			system(paste0("rm ",gsub(".tar.gz","",T_results[i])))
			
			system(paste0("tar zxvf ",AD_dir,"/2_Pvalue/",T_results[i]))
			tmp_AD_raw = read.table(gsub(".tar.gz","",T_results[i]),header=T)
			system(paste0("rm ",gsub(".tar.gz","",T_results[i])))
			
			thrgenome = -log10(0.05/nrow(tmp_T)); thrsuggest = -log10(1/nrow(tmp_T))
			
			if(max(tmp_T[,"MVN_logP_Minor"]) > thrsuggest){
				tmp_result_all2 = c()
				tmp_T_Sig = subset(tmp_T,tmp_T[,"MVN_logP_Minor"]  > thrsuggest)
				for(j in 1:length(unique(tmp_T_Sig[,"chr"]))){
					# [1] Get the T-statistics Model Info for the Peak SNP #
					tmp_T_Sig_chr = subset(tmp_T_Sig,tmp_T_Sig[,"chr"] == unique(tmp_T_Sig[,"chr"])[j])
					tmp_T_Sig_max = subset(tmp_T_Sig_chr,tmp_T_Sig_chr[,"MVN_logP_Minor"] == max(tmp_T_Sig_chr[,"MVN_logP_Minor"]))[1,]
					tmp_T_result = cbind(T_phes[i], tmp_T_Sig_max[,"chr"], tmp_T_Sig_max[,"pos"], tmp_T_Sig_max[,c("MVN_logP_Minor","tAB_AA","tAB_BB","betaAB_AA","betaAB_BB","varAB_AA","varAB_BB","df_res")])
				
					# [2] Get the A+D Model Info for the Peak SNP #
					tmp_AD_Sig_max = subset(tmp_AD_raw,tmp_AD_raw[,"chr"]==tmp_T_Sig_max[,"chr"] & tmp_AD_raw[,"pos"]==tmp_T_Sig_max[,"pos"])
					tmp_AD_result = tmp_AD_Sig_max[,c("X.logP.AddDomCode_Add.","X.logP.AddCode.","X.logP.DomCode.","X.logP.AddDomCode_Null.","X.logP.IndiCode.","betaAdd","betaDom")]
					tmp_AD_result = cbind(tmp_AD_result, tmp_AD_Sig_max[,"betaDom"]/tmp_AD_Sig_max[,"betaAdd"])
				
					# [3] Get the Confidence Level for the Peak SNP #
					tmp_Sig_SNP = subset(Freq_Info,Freq_Info[,1]==tmp_T_Sig_max[,"chr"] & Freq_Info[,4]==tmp_T_Sig_max[,"pos"])[,"SNP"]
					system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --r2 --ld-snp ",tmp_Sig_SNP," --ld-window 99999 --ld-window-kb 50000 --ld-window-r2 0"," --out tmp_Sig_SNP"))
					TEMP.ld <- read.table(file="tmp_Sig_SNP.ld",header=T,stringsAsFactors=F); rownames(TEMP.ld) = TEMP.ld[,"SNP_B"]; system("rm tmp_Sig_SNP.*")
					TEMP.ld.tmp = TEMP.ld[TEMP.ld[,"R2"]>CI_Pvalue_r2,]
					tmp_CI_result1 = paste(TEMP.ld.tmp[1,"BP_B"],TEMP.ld.tmp[nrow(TEMP.ld.tmp),"BP_B"],sep="-")
					
					tmp_T_chr = subset(tmp_T,tmp_T[,"chr"] == unique(tmp_T_Sig[,"chr"])[j] & tmp_T[,"MVN_logP_Minor"] >= (tmp_T_Sig_max[,"MVN_logP_Minor"]-CI_Pvalue_decrease))
					tmp_CI_result2 = paste(tmp_T_chr[1,"pos"],tmp_T_chr[nrow(tmp_T_chr),"pos"],sep="-")
				
					# [4] Get the MAF for the Peak SNP #
					tmp_MAF_result = subset(Freq_Info,Freq_Info[,1]==tmp_T_Sig_max[,"chr"] & Freq_Info[,4]==tmp_T_Sig_max[,"pos"])[,"MAF"]
					
					tmp_result = cbind(tmp_T_result,tmp_AD_result,tmp_MAF_result,tmp_CI_result1,tmp_CI_result2,thrgenome,thrsuggest)
					tmp_result_all2 = rbind(tmp_result_all2,tmp_result)
				}
				colnames(tmp_result_all2) = c("trait","chr","pos","MVN_logP_Minor","tAB_AA","tAB_BB","betaAB_AA","betaAB_BB","varAB_AA","varAB_BB","df_res","logP_AD_A","logP_A","logP_D","logP_AD","logP_Indi","beta_Add","beta_OverDom","OD/A","MAF","CI_r2","CI_logP","thrgenome","thrsuggest")
			
				print(paste0("Phe",i," : ",AD_phes[i]," has significant SNP in ",length(unique(tmp_T_Sig[,"chr"]))," chrs!"))
			}else{
				tmp_result_all2 = c()
				print(paste0("Phe",i," : ",T_phes[i]," no significant SNP!" ))
			}

			result_all2 = rbind(result_all2,tmp_result_all2)
		})
		}
		result_all2 = subset(result_all2,!(result_all2[,"trait"] %in% phe_rm) & !(result_all2[,"chr"] %in% chr_rm))
		write.table(result_all2,file="PeakSNP_T.txt",row.names=F,col.names=T,quote=F)

	}

	# STEP3: Filter overlap QTLs in two methods and SNPs with biggest P-value that are both >thrsuggest #

	result_all1 = read.table(paste0(AD_dir,"/PeakSNP_D.txt"),header=T)
	result_all2 = read.table(paste0(T_dir,"/PeakSNP_T.txt"),header=T)

	xxx = cbind(result_all2[,1:3],do.call(rbind,strsplit(result_all2[,"CI_r2"],"-"))); colnames(xxx) = c("trait","chr","pos","CI_START","CI_END") # Use xxx stands for T results #
	xxx[,"CI_START"] = as.numeric(xxx[,"CI_START"]); xxx[,"CI_END"] = as.numeric(xxx[,"CI_END"])
	yyy = cbind(result_all1[,1:3],do.call(rbind,strsplit(result_all1[,"CI_r2"],"-"))); colnames(yyy) = c("trait","chr","pos","CI_START","CI_END") # Use yyy stands for A+D results #
	yyy[,"CI_START"] = as.numeric(yyy[,"CI_START"]); yyy[,"CI_END"] = as.numeric(yyy[,"CI_END"])

	xxx = xxx[!duplicated(paste(xxx[,"trait"],xxx[,"chr"],xxx[,"pos"],sep=":")),]
	yyy = yyy[!duplicated(paste(yyy[,"trait"],yyy[,"chr"],yyy[,"pos"],sep=":")),]

	# [1] Get the overlap QTLs #
	tmp_result_all = data.frame()
	for(i in 1:nrow(yyy)){ # Based on A+D results, to search the overlap QTLs from T-statistics results #

		tmp_xxx = subset(xxx,xxx[,"trait"] == yyy[i,"trait"] & xxx[,"chr"] == yyy[i,"chr"])
		if(dim(tmp_xxx)[1] == 0){
			tmp_result = NA
		}else{	
			tmp_result = c()
			for(j in 1:nrow(tmp_xxx)){
				if( ((tmp_xxx[j,"CI_START"] >= yyy[i,"CI_START"]-Distance) & (tmp_xxx[j,"CI_START"] <= yyy[i,"CI_START"]+Distance)) | ((tmp_xxx[j,"CI_END"] >= yyy[i,"CI_END"]-Distance) & (tmp_xxx[j,"CI_END"] <= yyy[i,"CI_END"]+Distance)) ){
					tmp_res = paste(tmp_xxx[j,],collapse=":")
				}else{
					tmp_res = NA
				}
				tmp_result = c(tmp_result,tmp_res)
			}
		}
		
		tmp_result_all = rbind(tmp_result_all,tmp_result)
		print(i)

	}
	tmp_result_all = cbind(paste(yyy[,"trait"],yyy[,"chr"],yyy[,"pos"],yyy[,"CI_START"],yyy[,"CI_END"],sep=":"), tmp_result_all)
	colnames(tmp_result_all) = c("PeakSNP_D","PeakSNP_T")
	result_all = tmp_result_all[complete.cases(tmp_result_all),]

	# [2] Get the common significant SNPs (just choose the peak one) for each pair #
	result_all = cbind(result_all,Common_Peak=NA); result_detail = c()
	for(i in 1:nrow(result_all)){

		tmp_qtl1 = strsplit(result_all[i,1],":")[[1]]
		tmp_qtl2 = strsplit(result_all[i,2],":")[[1]]
		names(tmp_qtl1) = names(tmp_qtl2) = c("trait","chr","pos","CI_START","CI_END")
		
		setwd(AD_dir); system(paste0("tar zxvf 2_Pvalue/Pvalue_",tmp_qtl1[1],".txt.tar.gz"))
		tmp_pvalue1 = read.table(paste0("Pvalue_",tmp_qtl1[1],".txt"),header=T)
		tmp_pvalue1 = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]!=0 & tmp_pvalue1[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_pvalue1 = tmp_pvalue1[complete.cases(tmp_pvalue1),] # Remove SNPs with NA results #
		thrgenome1 = -log10(0.05/nrow(tmp_pvalue1)); thrsuggest1 = -log10(1/nrow(tmp_pvalue1))
		system(paste0("rm Pvalue_",tmp_qtl1[1],".txt"))
		
		setwd(T_dir); system(paste0("tar zxvf 2_Pvalue/Pvalue_",tmp_qtl2[1],".txt.tar.gz"))
		tmp_pvalue2 = read.table(paste0("Pvalue_",tmp_qtl2[1],".txt"),header=T)
		tmp_pvalue2 = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]!=0 & tmp_pvalue2[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_pvalue2 = tmp_pvalue2[complete.cases(tmp_pvalue2),] # Remove SNPs with NA results #
		thrgenome2 = -log10(0.05/nrow(tmp_pvalue2)); thrsuggest2 = -log10(1/nrow(tmp_pvalue2))
		system(paste0("rm Pvalue_",tmp_qtl2[1],".txt"))	
		
		if(tmp_qtl1[3] == tmp_qtl2[3]){ # Case1: If tmp_qtl1 and tmp_qtl2 share the common peak SNP #
			tmp_qtl = tmp_qtl1; result_all[i,3] = result_all[i,1]
			result_detail_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue1[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
			result_detail_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue2[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
			result_detail = rbind(result_detail,cbind(trait=tmp_qtl1[1],result_detail_AD,result_detail_T[-(1:2)]))
		}else{ # Case2: If tmp_qtl1 and tmp_qtl2 just have overlap CI region #
			
			region_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl1["chr"]) & tmp_pvalue1[,"pos"]>=as.numeric(tmp_qtl1["CI_START"]) & tmp_pvalue1[,"pos"]<=as.numeric(tmp_qtl1["CI_END"])) # Get the region SNPs for tmp_qtl1 #
			region_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl2["chr"]) & tmp_pvalue2[,"pos"]>=as.numeric(tmp_qtl2["CI_START"]) & tmp_pvalue2[,"pos"]<=as.numeric(tmp_qtl2["CI_END"]))  # Get the region SNPs for tmp_qtl2 #
			
			common_SNP = region_AD[,"pos"][region_AD[,"pos"] %in% region_T[,"pos"]]
			region_AD_common = region_AD[region_AD[,"pos"] %in% common_SNP,]
			region_AD_common_Sig = region_AD_common[region_AD_common[,"X.logP.DomCode."]>thrsuggest1,] # Get the common region & significant SNPs for tmp_qtl1 #
			region_T_common = region_T[region_T[,"pos"] %in% common_SNP,]
			region_T_common_Sig = region_T_common[region_T_common[,"MVN_logP_Minor"]>thrsuggest2,] # Get the common region & significant SNPs for tmp_qtl2 #
			
			if(nrow(region_AD_common_Sig) != 0 && nrow(region_T_common_Sig) != 0){
				common_SNP_Sig = region_AD_common_Sig[,"pos"][region_AD_common_Sig[,"pos"] %in% region_T_common_Sig[,"pos"]]
				if(length(common_SNP_Sig) > 0){
					region_AD_common_Key = region_AD_common_Sig[region_AD_common_Sig[,"pos"] %in% common_SNP_Sig,]
					tmp_qtl =c(tmp_qtl1[1], as.character(region_AD_common_Key[region_AD_common_Key[,"X.logP.DomCode."]==max(region_AD_common_Key[,"X.logP.DomCode."]),1:2][1,]))
					names(tmp_qtl) = c("trait","chr","pos")
					result_all[i,3] = paste(tmp_qtl,collapse=":")
					result_detail_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue1[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
					result_detail_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue2[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
					result_detail = rbind(result_detail,cbind(trait=tmp_qtl1[1],result_detail_AD,result_detail_T[-(1:2)]))
				}
			}

		}
		print(i)
	}	

	setwd(paste0(AD_dir,"/../Overlap-AD_T"))
	write.table(result_all,file="PeakSNP-D_T-Common1.txt",row.names=F,col.names=T,quote=F) # Contains all the overlap QTL-eQTL pairs and the 3rd col is the "Both-Sig SNP" #
	write.table(result_detail,file="PeakSNP-D_T-Common2.txt",row.names=F,col.names=T,quote=F) # Just contains the "Both-Sig SNPs" from overlap QTL-eQTL pairs #

})

##################################################################################################
### Part2: PeakSNP-AvsAD_T.r (Based on AvsAD P-value to choose the QTLs from A+D Model)
##################################################################################################

try({

	setwd(AD_dir)
	system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --recode --out tmp"))
	system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --freq --out tmp_freq"))
	Map_Info = read.table("tmp.map",header=F); Freq_Info = read.table("tmp_freq.frq",header=T)
	Freq_Info = cbind(Map_Info,Freq_Info)
	system("rm tmp*")

	AD_results = dir(paste0(AD_dir,"/2_Pvalue"),pattern=".txt.tar.gz")
	AD_phes = gsub("Pvalue_","",gsub(".txt.tar.gz","",AD_results)) 
	T_results = dir(paste0(T_dir,"/2_Pvalue"),pattern=".txt.tar.gz")
	T_phes = gsub("Pvalue_","",gsub(".txt.tar.gz","",T_results)) 

	# STEP1: Extract potential QTLs>thrsuggest based on A+D Results #
	setwd(AD_dir); result_all1 = c()
	for(i in 1:length(AD_phes)){
	try({

		system(paste0("tar zxvf 2_Pvalue/",AD_results[i]))
		tmp_AD_raw = read.table(gsub(".tar.gz","",AD_results[i]),header=T)
		tmp_AD = subset(tmp_AD_raw,tmp_AD_raw[,"chr"]!=0 & tmp_AD_raw[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_AD = tmp_AD[complete.cases(tmp_AD),] # Remove SNPs with NA results #
		system(paste0("rm ",gsub(".tar.gz","",AD_results[i])))
		
		system(paste0("tar zxvf ",T_dir,"/2_Pvalue/",AD_results[i]))
		tmp_T_raw = read.table(gsub(".tar.gz","",AD_results[i]),header=T)
		system(paste0("rm ",gsub(".tar.gz","",AD_results[i])))
		
		thrgenome = -log10(0.05/nrow(tmp_AD)); thrsuggest = -log10(1/nrow(tmp_AD))
		
		if(max(tmp_AD[,"X.logP.AddDomCode_Add."]) > thrsuggest){
			tmp_result_all1 = c()
			tmp_AD_Sig = subset(tmp_AD,tmp_AD[,"X.logP.AddDomCode_Add."]  > thrsuggest)
			for(j in 1:length(unique(tmp_AD_Sig[,"chr"]))){
				# [1] Get the A+D Model Info for the Peak SNP #
				tmp_AD_Sig_chr = subset(tmp_AD_Sig,tmp_AD_Sig[,"chr"] == unique(tmp_AD_Sig[,"chr"])[j])
				tmp_AD_Sig_max = subset(tmp_AD_Sig_chr,tmp_AD_Sig_chr[,"X.logP.AddDomCode_Add."] == max(tmp_AD_Sig_chr[,"X.logP.AddDomCode_Add."]))[1,]
				tmp_AD_result = cbind(AD_phes[i], tmp_AD_Sig_max[,"chr"], tmp_AD_Sig_max[,"pos"], tmp_AD_Sig_max[,c("X.logP.AddDomCode_Add.","X.logP.AddCode.","X.logP.DomCode.","X.logP.AddDomCode_Null.","X.logP.IndiCode.","betaAdd","betaDom")], tmp_AD_Sig_max[,"betaDom"]/tmp_AD_Sig_max[,"betaAdd"])
			
				# [2] Get the T-statistics Model Info for the Peak SNP #
				tmp_T_Sig_max = subset(tmp_T_raw,tmp_T_raw[,"chr"]==tmp_AD_Sig_max[,"chr"] & tmp_T_raw[,"pos"]==tmp_AD_Sig_max[,"pos"])
				tmp_T_result = tmp_T_Sig_max[,c("MVN_logP_Minor","tAB_AA","tAB_BB","betaAB_AA","betaAB_BB","varAB_AA","varAB_BB","df_res")]
			
				# [3] Get the Confidence Level for the Peak SNP #
				tmp_Sig_SNP = subset(Freq_Info,Freq_Info[,1]==tmp_AD_Sig_max[,"chr"] & Freq_Info[,4]==tmp_AD_Sig_max[,"pos"])[,"SNP"]
				system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --r2 --ld-snp ",tmp_Sig_SNP," --ld-window 99999 --ld-window-kb 50000 --ld-window-r2 0"," --out tmp_Sig_SNP"))
				TEMP.ld <- read.table(file="tmp_Sig_SNP.ld",header=T,stringsAsFactors=F); rownames(TEMP.ld) = TEMP.ld[,"SNP_B"]; system("rm tmp_Sig_SNP.*")
				TEMP.ld.tmp = TEMP.ld[TEMP.ld[,"R2"]>CI_Pvalue_r2,]
				tmp_CI_result1 = paste(TEMP.ld.tmp[1,"BP_B"],TEMP.ld.tmp[nrow(TEMP.ld.tmp),"BP_B"],sep="-")
				
				tmp_AD_chr = subset(tmp_AD,tmp_AD[,"chr"] == unique(tmp_AD_Sig[,"chr"])[j] & tmp_AD[,"X.logP.AddDomCode_Add."] >= (tmp_AD_Sig_max[,"X.logP.AddDomCode_Add."]-CI_Pvalue_decrease))
				tmp_CI_result2 = paste(tmp_AD_chr[1,"pos"],tmp_AD_chr[nrow(tmp_AD_chr),"pos"],sep="-")
			
				# [4] Get the MAF for the Peak SNP #
				tmp_MAF_result = subset(Freq_Info,Freq_Info[,1]==tmp_AD_Sig_max[,"chr"] & Freq_Info[,4]==tmp_AD_Sig_max[,"pos"])[,"MAF"]
				
				tmp_result = cbind(tmp_AD_result,tmp_T_result,tmp_MAF_result,tmp_CI_result1,tmp_CI_result2,thrgenome,thrsuggest)
				tmp_result_all1 = rbind(tmp_result_all1,tmp_result)
			}
			colnames(tmp_result_all1) = c("trait","chr","pos","logP_AD_A","logP_A","logP_D","logP_AD","logP_Indi","beta_Add","beta_OverDom","OD/A","MVN_logP_Minor","tAB_AA","tAB_BB","betaAB_AA","betaAB_BB","varAB_AA","varAB_BB","df_res","MAF","CI_r2","CI_logP","thrgenome","thrsuggest")
		
			print(paste0("Phe",i," : ",AD_phes[i]," has significant SNP in ",length(unique(tmp_AD_Sig[,"chr"]))," chrs!"))
		}else{
			tmp_result_all1 = c()
			print(paste0("Phe",i," : ",AD_phes[i]," no significant SNP." ))
		}
		
		result_all1 = rbind(result_all1,tmp_result_all1)
		
	})
	}
	result_all1 = subset(result_all1,!(result_all1[,"trait"] %in% phe_rm) & !(result_all1[,"chr"] %in% chr_rm))
	write.table(result_all1,file="PeakSNP_AvsAD.txt",row.names=F,col.names=T,quote=F)

	# STEP3: Filter overlap QTLs in two methods and SNPs with biggest P-value that are both >thrsuggest #

	result_all1 = read.table(paste0(AD_dir,"/PeakSNP_AvsAD.txt"),header=T)
	result_all2 = read.table(paste0(T_dir,"/PeakSNP_T.txt"),header=T)

	xxx = cbind(result_all2[,1:3],do.call(rbind,strsplit(result_all2[,"CI_r2"],"-"))); colnames(xxx) = c("trait","chr","pos","CI_START","CI_END") # Use xxx stands for T results #
	xxx[,"CI_START"] = as.numeric(xxx[,"CI_START"]); xxx[,"CI_END"] = as.numeric(xxx[,"CI_END"])
	yyy = cbind(result_all1[,1:3],do.call(rbind,strsplit(result_all1[,"CI_r2"],"-"))); colnames(yyy) = c("trait","chr","pos","CI_START","CI_END") # Use yyy stands for A+D results #
	yyy[,"CI_START"] = as.numeric(yyy[,"CI_START"]); yyy[,"CI_END"] = as.numeric(yyy[,"CI_END"])

	xxx = xxx[!duplicated(paste(xxx[,"trait"],xxx[,"chr"],xxx[,"pos"],sep=":")),]
	yyy = yyy[!duplicated(paste(yyy[,"trait"],yyy[,"chr"],yyy[,"pos"],sep=":")),]

	# [1] Get the overlap QTLs #
	tmp_result_all = data.frame()
	for(i in 1:nrow(yyy)){ # Based on A+D results, to search the overlap QTLs from T-statistics results #

		tmp_xxx = subset(xxx,xxx[,"trait"] == yyy[i,"trait"] & xxx[,"chr"] == yyy[i,"chr"])
		if(dim(tmp_xxx)[1] == 0){
			tmp_result = NA
		}else{	
			tmp_result = c()
			for(j in 1:nrow(tmp_xxx)){
				if( ((tmp_xxx[j,"CI_START"] >= yyy[i,"CI_START"]-Distance) & (tmp_xxx[j,"CI_START"] <= yyy[i,"CI_START"]+Distance)) | ((tmp_xxx[j,"CI_END"] >= yyy[i,"CI_END"]-Distance) & (tmp_xxx[j,"CI_END"] <= yyy[i,"CI_END"]+Distance)) ){
					tmp_res = paste(tmp_xxx[j,],collapse=":")
				}else{
					tmp_res = NA
				}
				tmp_result = c(tmp_result,tmp_res)
			}
		}
		
		tmp_result_all = rbind(tmp_result_all,tmp_result)
		print(i)

	}
	tmp_result_all = cbind(paste(yyy[,"trait"],yyy[,"chr"],yyy[,"pos"],yyy[,"CI_START"],yyy[,"CI_END"],sep=":"), tmp_result_all)
	colnames(tmp_result_all) = c("PeakSNP_AvsAD","PeakSNP_T")
	result_all = tmp_result_all[complete.cases(tmp_result_all),]

	# [2] Get the common significant SNPs (just choose the peak one) for each pair #
	result_all = cbind(result_all,Common_Peak=NA); result_detail = c()
	for(i in 1:nrow(result_all)){

		tmp_qtl1 = strsplit(result_all[i,1],":")[[1]]
		tmp_qtl2 = strsplit(result_all[i,2],":")[[1]]
		names(tmp_qtl1) = names(tmp_qtl2) = c("trait","chr","pos","CI_START","CI_END")
		
		setwd(AD_dir); system(paste0("tar zxvf 2_Pvalue/Pvalue_",tmp_qtl1[1],".txt.tar.gz"))
		tmp_pvalue1 = read.table(paste0("Pvalue_",tmp_qtl1[1],".txt"),header=T)
		tmp_pvalue1 = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]!=0 & tmp_pvalue1[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_pvalue1 = tmp_pvalue1[complete.cases(tmp_pvalue1),] # Remove SNPs with NA results #
		thrgenome1 = -log10(0.05/nrow(tmp_pvalue1)); thrsuggest1 = -log10(1/nrow(tmp_pvalue1))
		system(paste0("rm Pvalue_",tmp_qtl1[1],".txt"))
		
		setwd(T_dir); system(paste0("tar zxvf 2_Pvalue/Pvalue_",tmp_qtl2[1],".txt.tar.gz"))
		tmp_pvalue2 = read.table(paste0("Pvalue_",tmp_qtl2[1],".txt"),header=T)
		tmp_pvalue2 = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]!=0 & tmp_pvalue2[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_pvalue2 = tmp_pvalue2[complete.cases(tmp_pvalue2),] # Remove SNPs with NA results #
		thrgenome2 = -log10(0.05/nrow(tmp_pvalue2)); thrsuggest2 = -log10(1/nrow(tmp_pvalue2))
		system(paste0("rm Pvalue_",tmp_qtl2[1],".txt"))	
		
		if(tmp_qtl1[3] == tmp_qtl2[3]){ # Case1: If tmp_qtl1 and tmp_qtl2 share the common peak SNP #
			tmp_qtl = tmp_qtl1; result_all[i,3] = result_all[i,1]
			result_detail_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue1[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
			result_detail_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue2[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
			result_detail = rbind(result_detail,cbind(trait=tmp_qtl1[1],result_detail_AD,result_detail_T[-(1:2)]))
		}else{ # Case2: If tmp_qtl1 and tmp_qtl2 just have overlap CI region #
			
			region_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl1["chr"]) & tmp_pvalue1[,"pos"]>=as.numeric(tmp_qtl1["CI_START"]) & tmp_pvalue1[,"pos"]<=as.numeric(tmp_qtl1["CI_END"])) # Get the region SNPs for tmp_qtl1 #
			region_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl2["chr"]) & tmp_pvalue2[,"pos"]>=as.numeric(tmp_qtl2["CI_START"]) & tmp_pvalue2[,"pos"]<=as.numeric(tmp_qtl2["CI_END"]))  # Get the region SNPs for tmp_qtl2 #
			
			common_SNP = region_AD[,"pos"][region_AD[,"pos"] %in% region_T[,"pos"]]
			region_AD_common = region_AD[region_AD[,"pos"] %in% common_SNP,]
			region_AD_common_Sig = region_AD_common[region_AD_common[,"X.logP.AddDomCode_Add."]>thrsuggest1,] # Get the common region & significant SNPs for tmp_qtl1 #
			region_T_common = region_T[region_T[,"pos"] %in% common_SNP,]
			region_T_common_Sig = region_T_common[region_T_common[,"MVN_logP_Minor"]>thrsuggest2,] # Get the common region & significant SNPs for tmp_qtl2 #
			
			if(nrow(region_AD_common_Sig) != 0 && nrow(region_T_common_Sig) != 0){
				common_SNP_Sig = region_AD_common_Sig[,"pos"][region_AD_common_Sig[,"pos"] %in% region_T_common_Sig[,"pos"]]
				if(length(common_SNP_Sig) > 0){
					region_AD_common_Key = region_AD_common_Sig[region_AD_common_Sig[,"pos"] %in% common_SNP_Sig,]
					tmp_qtl =c(tmp_qtl1[1], as.character(region_AD_common_Key[region_AD_common_Key[,"X.logP.AddDomCode_Add."]==max(region_AD_common_Key[,"X.logP.AddDomCode_Add."]),1:2][1,]))
					names(tmp_qtl) = c("trait","chr","pos")
					result_all[i,3] = paste(tmp_qtl,collapse=":")
					result_detail_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue1[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
					result_detail_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue2[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
					result_detail = rbind(result_detail,cbind(trait=tmp_qtl1[1],result_detail_AD,result_detail_T[-(1:2)]))
				}
			}

		}
		print(i)
	}	

	setwd(paste0(AD_dir,"/../Overlap-AD_T"))
	write.table(result_all,file="PeakSNP-AvsAD_T-Common1.txt",row.names=F,col.names=T,quote=F) # Contains all the overlap QTL-eQTL pairs and the 3rd col is the "Both-Sig SNP" #
	write.table(result_detail,file="PeakSNP-AvsAD_T-Common2.txt",row.names=F,col.names=T,quote=F) # Just contains the "Both-Sig SNPs" from overlap QTL-eQTL pairs #


})

##################################################################################################
### Part3: PeakSNP-AD_T.r (Based on AD P-value to choose the QTLs from A+D Model)
##################################################################################################

try({
	
	setwd(AD_dir)
	system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --recode --out tmp"))
	system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --freq --out tmp_freq"))
	Map_Info = read.table("tmp.map",header=F); Freq_Info = read.table("tmp_freq.frq",header=T)
	Freq_Info = cbind(Map_Info,Freq_Info)
	system("rm tmp*")

	AD_results = dir(paste0(AD_dir,"/2_Pvalue"),pattern=".txt.tar.gz")
	AD_phes = gsub("Pvalue_","",gsub(".txt.tar.gz","",AD_results)) 
	T_results = dir(paste0(T_dir,"/2_Pvalue"),pattern=".txt.tar.gz")
	T_phes = gsub("Pvalue_","",gsub(".txt.tar.gz","",T_results)) 

	# STEP1: Extract potential QTLs>thrsuggest based on A+D Results #
	setwd(AD_dir); result_all1 = c()
	for(i in 1:length(AD_phes)){
	try({

		system(paste0("tar zxvf 2_Pvalue/",AD_results[i]))
		tmp_AD_raw = read.table(gsub(".tar.gz","",AD_results[i]),header=T)
		tmp_AD = subset(tmp_AD_raw,tmp_AD_raw[,"chr"]!=0 & tmp_AD_raw[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_AD = tmp_AD[complete.cases(tmp_AD),] # Remove SNPs with NA results #
		system(paste0("rm ",gsub(".tar.gz","",AD_results[i])))
		
		system(paste0("tar zxvf ",T_dir,"/2_Pvalue/",AD_results[i]))
		tmp_T_raw = read.table(gsub(".tar.gz","",AD_results[i]),header=T)
		system(paste0("rm ",gsub(".tar.gz","",AD_results[i])))
		
		thrgenome = -log10(0.05/nrow(tmp_AD)); thrsuggest = -log10(1/nrow(tmp_AD))
		
		if(max(tmp_AD[,"X.logP.AddDomCode_Null."]) > thrsuggest){
			tmp_result_all1 = c()
			tmp_AD_Sig = subset(tmp_AD,tmp_AD[,"X.logP.AddDomCode_Null."]  > thrsuggest)
			for(j in 1:length(unique(tmp_AD_Sig[,"chr"]))){
				# [1] Get the A+D Model Info for the Peak SNP #
				tmp_AD_Sig_chr = subset(tmp_AD_Sig,tmp_AD_Sig[,"chr"] == unique(tmp_AD_Sig[,"chr"])[j])
				tmp_AD_Sig_max = subset(tmp_AD_Sig_chr,tmp_AD_Sig_chr[,"X.logP.AddDomCode_Null."] == max(tmp_AD_Sig_chr[,"X.logP.AddDomCode_Null."]))[1,]
				tmp_AD_result = cbind(AD_phes[i], tmp_AD_Sig_max[,"chr"], tmp_AD_Sig_max[,"pos"], tmp_AD_Sig_max[,c("X.logP.AddDomCode_Add.","X.logP.AddCode.","X.logP.DomCode.","X.logP.AddDomCode_Null.","X.logP.IndiCode.","betaAdd","betaDom")], tmp_AD_Sig_max[,"betaDom"]/tmp_AD_Sig_max[,"betaAdd"])
			
				# [2] Get the T-statistics Model Info for the Peak SNP #
				tmp_T_Sig_max = subset(tmp_T_raw,tmp_T_raw[,"chr"]==tmp_AD_Sig_max[,"chr"] & tmp_T_raw[,"pos"]==tmp_AD_Sig_max[,"pos"])
				tmp_T_result = tmp_T_Sig_max[,c("MVN_logP_Minor","tAB_AA","tAB_BB","betaAB_AA","betaAB_BB","varAB_AA","varAB_BB","df_res")]
			
				# [3] Get the Confidence Level for the Peak SNP #
				tmp_Sig_SNP = subset(Freq_Info,Freq_Info[,1]==tmp_AD_Sig_max[,"chr"] & Freq_Info[,4]==tmp_AD_Sig_max[,"pos"])[,"SNP"]
				system(paste0("plink2 --noweb --tfile /SAN/mottlab/heterosis/",folder_name,"/1data/",file_name," --r2 --ld-snp ",tmp_Sig_SNP," --ld-window 99999 --ld-window-kb 50000 --ld-window-r2 0"," --out tmp_Sig_SNP"))
				TEMP.ld <- read.table(file="tmp_Sig_SNP.ld",header=T,stringsAsFactors=F); rownames(TEMP.ld) = TEMP.ld[,"SNP_B"]; system("rm tmp_Sig_SNP.*")
				TEMP.ld.tmp = TEMP.ld[TEMP.ld[,"R2"]>CI_Pvalue_r2,]
				tmp_CI_result1 = paste(TEMP.ld.tmp[1,"BP_B"],TEMP.ld.tmp[nrow(TEMP.ld.tmp),"BP_B"],sep="-")
				
				tmp_AD_chr = subset(tmp_AD,tmp_AD[,"chr"] == unique(tmp_AD_Sig[,"chr"])[j] & tmp_AD[,"X.logP.AddDomCode_Null."] >= (tmp_AD_Sig_max[,"X.logP.AddDomCode_Null."]-CI_Pvalue_decrease))
				tmp_CI_result2 = paste(tmp_AD_chr[1,"pos"],tmp_AD_chr[nrow(tmp_AD_chr),"pos"],sep="-")
			
				# [4] Get the MAF for the Peak SNP #
				tmp_MAF_result = subset(Freq_Info,Freq_Info[,1]==tmp_AD_Sig_max[,"chr"] & Freq_Info[,4]==tmp_AD_Sig_max[,"pos"])[,"MAF"]
				
				tmp_result = cbind(tmp_AD_result,tmp_T_result,tmp_MAF_result,tmp_CI_result1,tmp_CI_result2,thrgenome,thrsuggest)
				tmp_result_all1 = rbind(tmp_result_all1,tmp_result)
			}
			colnames(tmp_result_all1) = c("trait","chr","pos","logP_AD_A","logP_A","logP_D","logP_AD","logP_Indi","beta_Add","beta_OverDom","OD/A","MVN_logP_Minor","tAB_AA","tAB_BB","betaAB_AA","betaAB_BB","varAB_AA","varAB_BB","df_res","MAF","CI_r2","CI_logP","thrgenome","thrsuggest")
		
			print(paste0("Phe",i," : ",AD_phes[i]," has significant SNP in ",length(unique(tmp_AD_Sig[,"chr"]))," chrs!"))
		}else{
			tmp_result_all1 = c()
			print(paste0("Phe",i," : ",AD_phes[i]," no significant SNP." ))
		}
		
		result_all1 = rbind(result_all1,tmp_result_all1)
		
	})
	}
	result_all1 = subset(result_all1,!(result_all1[,"trait"] %in% phe_rm) & !(result_all1[,"chr"] %in% chr_rm))
	write.table(result_all1,file="PeakSNP_AD.txt",row.names=F,col.names=T,quote=F)

	# STEP3: Filter overlap QTLs in two methods and SNPs with biggest P-value that are both >thrsuggest #

	result_all1 = read.table(paste0(AD_dir,"/PeakSNP_AD.txt"),header=T)
	result_all2 = read.table(paste0(T_dir,"/PeakSNP_T.txt"),header=T)

	xxx = cbind(result_all2[,1:3],do.call(rbind,strsplit(result_all2[,"CI_r2"],"-"))); colnames(xxx) = c("trait","chr","pos","CI_START","CI_END") # Use xxx stands for T results #
	xxx[,"CI_START"] = as.numeric(xxx[,"CI_START"]); xxx[,"CI_END"] = as.numeric(xxx[,"CI_END"])
	yyy = cbind(result_all1[,1:3],do.call(rbind,strsplit(result_all1[,"CI_r2"],"-"))); colnames(yyy) = c("trait","chr","pos","CI_START","CI_END") # Use yyy stands for A+D results #
	yyy[,"CI_START"] = as.numeric(yyy[,"CI_START"]); yyy[,"CI_END"] = as.numeric(yyy[,"CI_END"])

	xxx = xxx[!duplicated(paste(xxx[,"trait"],xxx[,"chr"],xxx[,"pos"],sep=":")),]
	yyy = yyy[!duplicated(paste(yyy[,"trait"],yyy[,"chr"],yyy[,"pos"],sep=":")),]

	# [1] Get the overlap QTLs #
	tmp_result_all = data.frame()
	for(i in 1:nrow(yyy)){ # Based on A+D results, to search the overlap QTLs from T-statistics results #

		tmp_xxx = subset(xxx,xxx[,"trait"] == yyy[i,"trait"] & xxx[,"chr"] == yyy[i,"chr"])
		if(dim(tmp_xxx)[1] == 0){
			tmp_result = NA
		}else{	
			tmp_result = c()
			for(j in 1:nrow(tmp_xxx)){
				if( ((tmp_xxx[j,"CI_START"] >= yyy[i,"CI_START"]-Distance) & (tmp_xxx[j,"CI_START"] <= yyy[i,"CI_START"]+Distance)) | ((tmp_xxx[j,"CI_END"] >= yyy[i,"CI_END"]-Distance) & (tmp_xxx[j,"CI_END"] <= yyy[i,"CI_END"]+Distance)) ){
					tmp_res = paste(tmp_xxx[j,],collapse=":")
				}else{
					tmp_res = NA
				}
				tmp_result = c(tmp_result,tmp_res)
			}
		}
		
		tmp_result_all = rbind(tmp_result_all,tmp_result)
		print(i)

	}
	tmp_result_all = cbind(paste(yyy[,"trait"],yyy[,"chr"],yyy[,"pos"],yyy[,"CI_START"],yyy[,"CI_END"],sep=":"), tmp_result_all)
	colnames(tmp_result_all) = c("PeakSNP_AD","PeakSNP_T")
	result_all = tmp_result_all[complete.cases(tmp_result_all),]

	# [2] Get the common significant SNPs (just choose the peak one) for each pair #
	result_all = cbind(result_all,Common_Peak=NA); result_detail = c()
	for(i in 1:nrow(result_all)){

		tmp_qtl1 = strsplit(result_all[i,1],":")[[1]]
		tmp_qtl2 = strsplit(result_all[i,2],":")[[1]]
		names(tmp_qtl1) = names(tmp_qtl2) = c("trait","chr","pos","CI_START","CI_END")
		
		setwd(AD_dir); system(paste0("tar zxvf 2_Pvalue/Pvalue_",tmp_qtl1[1],".txt.tar.gz"))
		tmp_pvalue1 = read.table(paste0("Pvalue_",tmp_qtl1[1],".txt"),header=T)
		tmp_pvalue1 = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]!=0 & tmp_pvalue1[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_pvalue1 = tmp_pvalue1[complete.cases(tmp_pvalue1),] # Remove SNPs with NA results #
		thrgenome1 = -log10(0.05/nrow(tmp_pvalue1)); thrsuggest1 = -log10(1/nrow(tmp_pvalue1))
		system(paste0("rm Pvalue_",tmp_qtl1[1],".txt"))
		
		setwd(T_dir); system(paste0("tar zxvf 2_Pvalue/Pvalue_",tmp_qtl2[1],".txt.tar.gz"))
		tmp_pvalue2 = read.table(paste0("Pvalue_",tmp_qtl2[1],".txt"),header=T)
		tmp_pvalue2 = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]!=0 & tmp_pvalue2[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_pvalue2 = tmp_pvalue2[complete.cases(tmp_pvalue2),] # Remove SNPs with NA results #
		thrgenome2 = -log10(0.05/nrow(tmp_pvalue2)); thrsuggest2 = -log10(1/nrow(tmp_pvalue2))
		system(paste0("rm Pvalue_",tmp_qtl2[1],".txt"))	
		
		if(tmp_qtl1[3] == tmp_qtl2[3]){ # Case1: If tmp_qtl1 and tmp_qtl2 share the common peak SNP #
			tmp_qtl = tmp_qtl1; result_all[i,3] = result_all[i,1]
			result_detail_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue1[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
			result_detail_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue2[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
			result_detail = rbind(result_detail,cbind(trait=tmp_qtl1[1],result_detail_AD,result_detail_T[-(1:2)]))
		}else{ # Case2: If tmp_qtl1 and tmp_qtl2 just have overlap CI region #
			
			region_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl1["chr"]) & tmp_pvalue1[,"pos"]>=as.numeric(tmp_qtl1["CI_START"]) & tmp_pvalue1[,"pos"]<=as.numeric(tmp_qtl1["CI_END"])) # Get the region SNPs for tmp_qtl1 #
			region_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl2["chr"]) & tmp_pvalue2[,"pos"]>=as.numeric(tmp_qtl2["CI_START"]) & tmp_pvalue2[,"pos"]<=as.numeric(tmp_qtl2["CI_END"]))  # Get the region SNPs for tmp_qtl2 #
			
			common_SNP = region_AD[,"pos"][region_AD[,"pos"] %in% region_T[,"pos"]]
			region_AD_common = region_AD[region_AD[,"pos"] %in% common_SNP,]
			region_AD_common_Sig = region_AD_common[region_AD_common[,"X.logP.AddDomCode_Null."]>thrsuggest1,] # Get the common region & significant SNPs for tmp_qtl1 #
			region_T_common = region_T[region_T[,"pos"] %in% common_SNP,]
			region_T_common_Sig = region_T_common[region_T_common[,"MVN_logP_Minor"]>thrsuggest2,] # Get the common region & significant SNPs for tmp_qtl2 #
			
			if(nrow(region_AD_common_Sig) != 0 && nrow(region_T_common_Sig) != 0){
				common_SNP_Sig = region_AD_common_Sig[,"pos"][region_AD_common_Sig[,"pos"] %in% region_T_common_Sig[,"pos"]]
				if(length(common_SNP_Sig) > 0){
					region_AD_common_Key = region_AD_common_Sig[region_AD_common_Sig[,"pos"] %in% common_SNP_Sig,]
					tmp_qtl =c(tmp_qtl1[1], as.character(region_AD_common_Key[region_AD_common_Key[,"X.logP.AddDomCode_Null."]==max(region_AD_common_Key[,"X.logP.AddDomCode_Null."]),1:2][1,]))
					names(tmp_qtl) = c("trait","chr","pos")
					result_all[i,3] = paste(tmp_qtl,collapse=":")
					result_detail_AD = subset(tmp_pvalue1,tmp_pvalue1[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue1[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
					result_detail_T = subset(tmp_pvalue2,tmp_pvalue2[,"chr"]==as.numeric(tmp_qtl["chr"]) & tmp_pvalue2[,"pos"]==as.numeric(tmp_qtl["pos"]))[1,]
					result_detail = rbind(result_detail,cbind(trait=tmp_qtl1[1],result_detail_AD,result_detail_T[-(1:2)]))
				}
			}

		}
		print(i)
	}	

	setwd(paste0(AD_dir,"/../Overlap-AD_T"))
	write.table(result_all,file="PeakSNP-AD_T-Common1.txt",row.names=F,col.names=T,quote=F) # Contains all the overlap QTL-eQTL pairs and the 3rd col is the "Both-Sig SNP" #
	write.table(result_detail,file="PeakSNP-AD_T-Common2.txt",row.names=F,col.names=T,quote=F) # Just contains the "Both-Sig SNPs" from overlap QTL-eQTL pairs #


})


