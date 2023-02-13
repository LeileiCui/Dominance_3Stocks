
### STEP1. Generate the scaled residual GenABEL ###

rm(list=ls());options(stringsAsFactors=FALSE)
setwd("/SAN/mottlab/heterosis/4_HSmice/5VC/data")

# [1] Load the clean GenABEL data #
library(GenABEL)
load("Data_clean.ABEL.dat"); dat = dat_clean 
export.plink(dat,filebasename="Data_clean",phenotypes="all")

# [2] Calculate the residuals of each phenotype by control age & batch #
pheno_all = read.table(file="Data_clean.phe",header=T)[,-1]
rownames(pheno_all) = 1:nrow(pheno_all) # In order to have the same rowname with "fit$residuals" #

pheno = pheno_all[,-c(1:2)] # The first 11 cols are: IID and 1 covariates #

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

# [3] Scale each residual by mean & sd #
pheno_all[,-c(1:2)] = pheno # Revalue the 253 phenotypes by the corresponding residual values #
mean = apply(pheno_all[,-c(1:2)],2,mean,na.rm=T)
sd = apply(pheno_all[,-c(1:2)],2,sd,na.rm=T)
pheno_all[,-c(1:2)] = sweep(pheno_all[,-c(1:2)],2,mean,"-")
pheno_all[,-c(1:2)] = sweep(pheno_all[,-c(1:2)],2,sd,"/")

# [4] Revalue the phdata of GenABEL data and output #
phdata(dat) = pheno_all
save(dat,file="HSmice_res.ABEL.dat")
export.plink(dat,"HSmice_res",phenotypes=NULL)

HSmice.fam = read.table(file="HSmice_res.tfam",header=F,stringsAsFactors=F)
HSmice.phen = data.frame(HSmice.fam[,1:2],pheno)
write.table(HSmice.phen,file="HSmice_res.pheno",col.names=F,quote=F,row.names=F)
write.table(colnames(pheno),file="HSmice_res.traitnames",col.names=F,quote=F,row.names=F)

	# RemoveInd = read.table(file="HSmice_Info_RemoveInd.txt",header=F)
	# write.table(HSmice.phen[-c(RemoveInd[,1]),],file="HSmice_resrm.pheno",col.names=F,quote=F,row.names=F)
	# write.table(colnames(pheno),file="HSmice_resrm.traitnames",col.names=F,quote=F,row.names=F)

# RemovePhe = read.table(file="HSmice_Info_RemovePhe.txt",header=F)[,1]

# HSmice.fam = read.table(file="HSmice_res.tfam",header=F,stringsAsFactors=F)
# HSmice.phen = data.frame(HSmice.fam[,1:2],pheno)
# write.table(HSmice.phen[,!(colnames(HSmice.phen) %in% RemovePhe)],file="HSmice_res.pheno",col.names=F,quote=F,row.names=F)
# write.table(colnames(pheno)[!colnames(pheno) %in% RemovePhe],file="HSmice_res.traitnames",col.names=F,quote=F,row.names=F)

# RemoveInd = read.table(file="HSmice_Info_RemoveInd.txt",header=F)
# write.table(HSmice.phen[-c(RemoveInd[,1]),!(colnames(HSmice.phen) %in% RemovePhe)],file="HSmice_resrm.pheno",col.names=F,quote=F,row.names=F)
# write.table(colnames(pheno)[!colnames(pheno) %in% RemovePhe],file="HSmice_resrm.traitnames",col.names=F,quote=F,row.names=F)


