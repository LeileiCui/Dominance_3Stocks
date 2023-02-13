
### STEP1. Generate the scaled residual GenABEL ###

rm(list=ls());options(stringsAsFactors=FALSE)
setwd("/SAN/mottlab/heterosis/3_HSrats/5VC/data")

# [1] Load the clean GenABEL data #
library(GenABEL)
load("Data_clean.ABEL.dat")
#load('/SAN/mottlab/heterosis/3_HSrats/OLD-KEY/1data/2Oxford-Data/Key-Data/final_mat.RData')
load('/SAN/mottlab/heterosis/0_Data/2_HS-RAT/RESIDUALS/final_mat.RData')
final_mat[,1] = gsub(" ","",final_mat[,1]) # Convert "Batch 1" into "Batch1" #
#dat_Amelie = final_mat[rownames(dat_clean@phdata),gsub("_qn","",gsub("_bc","",colnames(dat_clean@phdata)[-1]))]
dat_Amelie = final_mat[rownames(dat_clean@phdata),colnames(dat_clean@phdata)[-1]]
phdata(dat_clean) = cbind("id"=rownames(dat_Amelie),dat_Amelie); dat = dat_clean
phenotype_Amelie = cbind("FID"=1:nrow(dat_Amelie),"IID"=rownames(dat_Amelie),dat_Amelie)
write.table(phenotype_Amelie,file="Data_clean.phe",row.names=F,col.names=T,quote=F)

# [2] Calculate the residuals of each phenotype by control age & batch #
pheno_all = read.table(file="Data_clean.phe",header=T)[,-1]
rownames(pheno_all) = 1:nrow(pheno_all) # In order to have the same rowname with "fit$residuals" #
pheno_all[,"batch"] = as.factor(pheno_all[,"batch"]) # In order to use "batch" as a random variable #

pheno = pheno_all[,-c(1:11)] # The first 11 cols are: IID and 10 covariates #
for(i in 1:ncol(pheno)){
    ageLevels = length(unique(pheno_all[!is.na(pheno[,i]),"age"]))
	batchLevels = length(unique(pheno_all[!is.na(pheno[,i]),"batch"]))
    sexLevels = length(unique(pheno_all[!is.na(pheno[,i]),"sex"]))
	
    if(ageLevels>1 & batchLevels>1 & sexLevels>1){
        fit = lm(pheno[,i]~pheno_all[,"age"]+pheno_all[,"batch"]+pheno_all[,"sex"]); idx = as.numeric(names(fit$residuals)); pheno[idx,i] = fit$residuals
    }
    if(ageLevels>1 & batchLevels>1 & sexLevels==1){
        fit = lm(pheno[,i]~pheno_all[,"age"]+pheno_all[,"batch"]); idx = as.numeric(names(fit$residuals)); pheno[idx,i] = fit$residuals
    }
	if(ageLevels>1 & batchLevels==1 & sexLevels>1){
        fit = lm(pheno[,i]~pheno_all[,"age"]+pheno_all[,"sex"]); idx = as.numeric(names(fit$residuals)); pheno[idx,i] = fit$residuals
    }
	if(ageLevels==1 & batchLevels>1 & sexLevels>1){
        fit = lm(pheno[,i]~pheno_all[,"batch"]+pheno_all[,"sex"]); idx = as.numeric(names(fit$residuals)); pheno[idx,i] = fit$residuals
    }
    if(ageLevels>1 & batchLevels==1 & sexLevels==1){
        fit = lm(pheno[,i]~pheno_all[,"age"]); idx = as.numeric(names(fit$residuals)); pheno[idx,i] = fit$residuals
    }
	if(ageLevels==1 & batchLevels>1 & sexLevels==1){
        fit = lm(pheno[,i]~pheno_all[,"batch"]); idx = as.numeric(names(fit$residuals)); pheno[idx,i] = fit$residuals
    }
	if(ageLevels==1 & batchLevels==1 & sexLevels>1){
        fit = lm(pheno[,i]~pheno_all[,"sex"]); idx = as.numeric(names(fit$residuals)); pheno[idx,i] = fit$residuals
    }
    cat(i,"\n")
}

# [3] Scale each residual by mean & sd #
pheno_all[,-c(1:11)] = pheno # Revalue the 253 phenotypes by the corresponding residual values #
mean = apply(pheno_all[,-c(1:11)],2,mean,na.rm=T)
sd = apply(pheno_all[,-c(1:11)],2,sd,na.rm=T)
pheno_all[,-c(1:11)] = sweep(pheno_all[,-c(1:11)],2,mean,"-")
pheno_all[,-c(1:11)] = sweep(pheno_all[,-c(1:11)],2,sd,"/")

# [4] Revalue the phdata of GenABEL data and output #
phdata(dat) = pheno_all
save(dat,file="HSrats_res.ABEL.dat")
export.plink(dat,"HSrats_res",phenotypes=NULL)

HSrats.fam = read.table(file="HSrats_res.tfam",header=F,stringsAsFactors=F)
HSrats.phen = data.frame(HSrats.fam[,1:2],pheno)
write.table(HSrats.phen,file="HSrats_res.pheno",col.names=F,quote=F,row.names=F)
write.table(colnames(pheno),file="HSrats_res.traitnames",col.names=F,quote=F,row.names=F)

	# RemoveInd = read.table(file="HSrats_Info_RemoveInd.txt",header=F)
	# write.table(HSrats.phen[-c(RemoveInd[,1]),],file="HSrats_resrm.pheno",col.names=F,quote=F,row.names=F)
	# write.table(colnames(pheno),file="HSrats_resrm.traitnames",col.names=F,quote=F,row.names=F)

# RemovePhe = read.table(file="HSrats_Info_RemovePhe.txt",header=F)[,1]

# HSrats.fam = read.table(file="HSrats_res.tfam",header=F,stringsAsFactors=F)
# HSrats.phen = data.frame(HSrats.fam[,1:2],pheno)
# write.table(HSrats.phen[,!(colnames(HSrats.phen) %in% RemovePhe)],file="HSrats_res.pheno",col.names=F,quote=F,row.names=F)
# write.table(colnames(pheno)[!colnames(pheno) %in% RemovePhe],file="HSrats_res.traitnames",col.names=F,quote=F,row.names=F)

# RemoveInd = read.table(file="HSrats_Info_RemoveInd.txt",header=F)
# write.table(HSrats.phen[-c(RemoveInd[,1]),!(colnames(HSrats.phen) %in% RemovePhe)],file="HSrats_resrm.pheno",col.names=F,quote=F,row.names=F)
# write.table(colnames(pheno)[!colnames(pheno) %in% RemovePhe],file="HSrats_resrm.traitnames",col.names=F,quote=F,row.names=F)

