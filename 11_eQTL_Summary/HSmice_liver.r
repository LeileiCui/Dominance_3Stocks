rm(list=ls());options(stringsAsFactors=FALSE)

out_dir = "/SAN/mottlab/heterosis/4_HSmice/3NonAdd_eQTL_Res/1_A+D_liver/SNP_Peak"
rawdat_dir = "/SAN/mottlab/heterosis/4_HSmice/1data/RAW_DGE/HSmice_liver"
pvalue_dir = "/SAN/mottlab/heterosis/4_HSmice/3NonAdd_eQTL_Res/1_A+D_liver/2_Pvalue"


for(pvalue_choose in 1:4){

	CI_Pvalue_decrease = 2; CI_Pvalue_r2 = 0.8
	pvalue_filter = c("X.logP.AddDomCode_Add.","X.logP.AddCode.","X.logP.DomCode.","X.logP.AddDomCode_Null.")[pvalue_choose]
	pvalue_outname = c("PeakSNP_AvsAD_thrsuggest.txt", "PeakSNP_A_thrsuggest.txt", "PeakSNP_D_thrsuggest.txt", "PeakSNP_AD_thrsuggest.txt")[pvalue_choose]

	setwd(out_dir)
	system(paste0("plink --bfile ", rawdat_dir, " --freq --out tmp_freq"))
	Map_Info = read.table(paste0(rawdat_dir, ".bim"),header=F); Freq_Info = read.table("tmp_freq.frq",header=T)
	Freq_Info = cbind(Map_Info[,1:4],Freq_Info)
	system("rm tmp*")

	AD_results = dir(pvalue_dir, pattern=".txt.tar.gz")
	AD_phes = gsub("Pvalue_","",gsub(".txt.tar.gz","",AD_results)) 

	# STEP1: Extract potential QTLs>thrsuggest based on A+D Results #
	setwd(out_dir); result_all1 = c()
	for(i in 1:length(AD_results)){
	try({

		system(paste0("tar zxvf ", pvalue_dir, "/", AD_results[i]))
		tmp_AD_raw = read.table(gsub(".tar.gz","",AD_results[i]),header=T)
		tmp_AD = subset(tmp_AD_raw,tmp_AD_raw[,"chr"]!=0 & tmp_AD_raw[,"pos"]!=0) # Remove SNPs without CHR and POS info #
		tmp_AD = tmp_AD[complete.cases(tmp_AD),] # Remove SNPs with NA results #
		system(paste0("rm ",gsub(".tar.gz","",AD_results[i])))
		
		thrgenome = -log10(0.05/nrow(tmp_AD)); thrsuggest = -log10(1/nrow(tmp_AD))
		
		if(max(tmp_AD[, pvalue_filter]) > thrsuggest){
			tmp_result_all1 = c()
			tmp_AD_Sig = subset(tmp_AD,tmp_AD[,pvalue_filter]  > thrsuggest)
			for(j in 1:length(unique(tmp_AD_Sig[,"chr"]))){
				# [1] Get the A+D Model Info for the Peak SNP #
				tmp_AD_Sig_chr = subset(tmp_AD_Sig,tmp_AD_Sig[,"chr"] == unique(tmp_AD_Sig[,"chr"])[j])
				tmp_AD_Sig_max = subset(tmp_AD_Sig_chr,tmp_AD_Sig_chr[,pvalue_filter] == max(tmp_AD_Sig_chr[,pvalue_filter]))[1,]
				tmp_AD_result = cbind(AD_phes[i], tmp_AD_Sig_max[,"chr"], tmp_AD_Sig_max[,"pos"], tmp_AD_Sig_max[,c("X.logP.AddDomCode_Add.","X.logP.AddCode.","X.logP.DomCode.","X.logP.AddDomCode_Null.","betaAdd","betaDom")], tmp_AD_Sig_max[,"betaDom"]/tmp_AD_Sig_max[,"betaAdd"])
			
				# [3] Get the Confidence Level for the Peak SNP #
				tmp_Sig_SNP = subset(Freq_Info,Freq_Info[,1]==tmp_AD_Sig_max[,"chr"] & Freq_Info[,4]==tmp_AD_Sig_max[,"pos"])[,"SNP"]
				system(paste0("plink --bfile ", rawdat_dir, " --r2 --ld-snp ",tmp_Sig_SNP," --ld-window 99999 --ld-window-kb 50000 --ld-window-r2 0"," --out tmp_Sig_SNP"))
				TEMP.ld <- read.table(file="tmp_Sig_SNP.ld",header=T,stringsAsFactors=F); rownames(TEMP.ld) = TEMP.ld[,"SNP_B"]; system("rm tmp_Sig_SNP.*")
				TEMP.ld.tmp = TEMP.ld[TEMP.ld[,"R2"]>CI_Pvalue_r2,]
				tmp_CI_result1 = paste(TEMP.ld.tmp[1,"BP_B"],TEMP.ld.tmp[nrow(TEMP.ld.tmp),"BP_B"],sep="-")
				
				tmp_AD_chr = subset(tmp_AD,tmp_AD[,"chr"] == unique(tmp_AD_Sig[,"chr"])[j] & tmp_AD[,pvalue_filter] >= (tmp_AD_Sig_max[,pvalue_filter]-CI_Pvalue_decrease))
				tmp_CI_result2 = paste(tmp_AD_chr[1,"pos"],tmp_AD_chr[nrow(tmp_AD_chr),"pos"],sep="-")
			
				# [4] Get the MAF for the Peak SNP #
				tmp_MAF_result = subset(Freq_Info,Freq_Info[,1]==tmp_AD_Sig_max[,"chr"] & Freq_Info[,4]==tmp_AD_Sig_max[,"pos"])[,"MAF"]
				
				tmp_result = cbind(tmp_AD_result,tmp_MAF_result,tmp_CI_result1,tmp_CI_result2,thrgenome,thrsuggest)
				tmp_result_all1 = rbind(tmp_result_all1,tmp_result)
			}
			colnames(tmp_result_all1) = c("trait","chr","pos","logP_AD_A","logP_A","logP_D","logP_AD","beta_Add","beta_OverDom","OD/A","MAF","CI_r2","CI_logP","thrgenome","thrsuggest")
		
			print(paste0("Phe",i," : ",AD_phes[i]," has significant SNP in ",length(unique(tmp_AD_Sig[,"chr"]))," chrs!"))
		}else{
			tmp_result_all1 = c()
			print(paste0("Phe",i," : ",AD_phes[i]," no significant SNP." ))
		}
		
		result_all1 = rbind(result_all1,tmp_result_all1)
		
	})
	}
	write.table(result_all1, file=pvalue_outname, row.names=F, col.names=T, quote=F)

}

