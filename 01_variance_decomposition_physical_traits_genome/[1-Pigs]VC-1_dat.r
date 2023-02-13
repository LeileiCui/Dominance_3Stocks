
### STEP1. Generate the scaled residual GenABEL ###

rm(list=ls());options(stringsAsFactors=FALSE)
setwd("/SAN/mottlab/heterosis/1_F2pigs/5VC/data")

# [1] Load the clean GenABEL data #
library(GenABEL)
load("Data_clean.ABEL.dat"); dat = dat_clean 
export.plink(dat,filebasename="Data_clean",phenotypes="all")

# [2] Calculate the residuals of each phenotype by control sex & batch #
pheno_all = read.table(file="Data_clean.phe",header=T)[,-1]
rownames(pheno_all) = 1:nrow(pheno_all) # In order to have the same rowname with "fit$residuals" #
pheno_all[,"batch"] = as.factor(pheno_all[,"batch"]) # In order to use "batch" as a random variable #

pheno = pheno_all[,-c(1:3)] # The first 3 cols are: IID;sex;batch #

for(i in 1:ncol(pheno)){
    batchLevels = length(unique(pheno_all[!is.na(pheno[,i]),"batch"]))
    sexLevels = length(unique(pheno_all[!is.na(pheno[,i]),"sex"]))
    if(batchLevels==1 & sexLevels>1){
        fit = lm(pheno[,i]~pheno_all[,"sex"])
        idx = as.numeric(names(fit$residuals))
        pheno[idx,i] = fit$residuals
    }
    else if(batchLevels>1 & sexLevels>1){
        fit = lm(pheno[,i]~pheno_all[,"sex"]+pheno_all[,"batch"])
        idx = as.numeric(names(fit$residuals))
        pheno[idx,i] = fit$residuals
    }
    else if(batchLevels>1 & sexLevels==1){
        fit = lm(pheno[,i]~pheno_all[,"batch"])
        idx = as.numeric(names(fit$residuals))
        pheno[idx,i] = fit$residuals
    }
    cat(i,"\n")
}

# [3] Scale each residual by mean & sd #

quantile_norm <- function(x){
try({
	names(x) = paste("Ind",1:length(x),sep="")
	x_clean = x[!is.na(x)]
	x_clean_norm = qnorm(rank(x_clean)/(1+length(x_clean)))
	x[names(x_clean)] = x_clean_norm
	return(x)
})
}
Normal_method = "QUANTILE"
for(i in 1:(dim(pheno)[2])){ if(Normal_method == "QUANTILE"){ pheno[,i] = quantile_norm(pheno[,i]) } }

pheno_all[,-c(1:3)] = pheno # Revalue the 253 phenotypes by the corresponding residual values #
mean = apply(pheno_all[,-c(1:3)],2,mean,na.rm=T)
sd = apply(pheno_all[,-c(1:3)],2,sd,na.rm=T)
pheno_all[,-c(1:3)] = sweep(pheno_all[,-c(1:3)],2,mean,"-")
pheno_all[,-c(1:3)] = sweep(pheno_all[,-c(1:3)],2,sd,"/")

# [4] Revalue the phdata of GenABEL data and output #
phdata(dat) = pheno_all
save(dat,file="F2_res.ABEL.dat")
export.plink(dat,"F2_res",phenotypes=NULL)

#F2.fam = read.table(file="F2_res.tfam",header=F,stringsAsFactors=F)
#F2.phen = data.frame(F2.fam[,1:2],pheno)
#write.table(F2.phen,file="F2_res.pheno",col.names=F,quote=F,row.names=F)
#write.table(colnames(pheno),file="F2_res.traitnames",col.names=F,quote=F,row.names=F)

#RemoveInd = read.table(file="F2_Info_RemoveInd.txt",header=F)
#write.table(F2.phen[-c(RemoveInd[,1]),],file="F2_resrm.pheno",col.names=F,quote=F,row.names=F)
#write.table(colnames(pheno),file="F2_resrm.traitnames",col.names=F,quote=F,row.names=F)

RemovePhe = read.table(file="F2_Info_RemovePhe.txt",header=F)[,1]

F2.fam = read.table(file="F2_res.tfam",header=F,stringsAsFactors=F)
F2.phen = data.frame(F2.fam[,1:2],pheno)
write.table(F2.phen[,!(colnames(F2.phen) %in% RemovePhe)],file="F2_res.pheno",col.names=F,quote=F,row.names=F)
write.table(colnames(pheno)[!colnames(pheno) %in% RemovePhe],file="F2_res.traitnames",col.names=F,quote=F,row.names=F)

RemoveInd = read.table(file="F2_Info_RemoveInd.txt",header=F)
write.table(F2.phen[-c(RemoveInd[,1]),!(colnames(F2.phen) %in% RemovePhe)],file="F2_resrm.pheno",col.names=F,quote=F,row.names=F)
write.table(colnames(pheno)[!colnames(pheno) %in% RemovePhe],file="F2_resrm.traitnames",col.names=F,quote=F,row.names=F)

