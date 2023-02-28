
#########################################################################################################
### STEP1. Generate "GT_Minimum_all.RData" including the Minimum info of all loci across all traits 
#########################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

library(parallel); library(data.table); num_nodes = 10
setwd("/SAN/mottlab/heterosis/4_HSmice/3NonAdd_QTL/1_A+D/1_PheGen/0_Data")

phe_all = read.table("/SAN/mottlab/heterosis/4_HSmice/3NonAdd_QTL/1_A+D/1_PheGen/5_PheTrans/Phe_Trans.txt",header=T)
read <- function(...) as.data.frame(fread(header=F,colClasses="character",...)); Chip_tfam <- read("Data_clean.tfam"); rownames(Chip_tfam) = Chip_tfam[,2]
Chip_map = read.table("Data_clean.map",header=F)

GT_sum_all = as.data.frame(matrix(NA,dim(Chip_map)[1],(dim(phe_all)[2]-3)))
rownames(GT_sum_all) = Chip_map[,2]; colnames(GT_sum_all) = colnames(phe_all)[-(1:3)]
for(i in 3:length(colnames(phe_all))){
try({
	
	phenotype_name = colnames(phe_all)[i]
	phe_tmp = phe_all[,c(1,i)][complete.cases(phe_all[,c(1,i)]),]; keep_inds = rownames(phe_tmp)
	
	myfolder=phenotype_name; if(myfolder%in%dir()==FALSE) dir.create(myfolder)
	write.table(Chip_tfam[keep_inds,1:2],file=paste(myfolder,"/TMP_indlist.txt",sep=""),row.names=F,col.names=F,quote=F)
	system(paste("plink2 --noweb --tfile Data_clean --keep ",myfolder,"/TMP_indlist.txt --recode --out ",myfolder,"/TMP_indlist",sep=""))
	Chip_map <- read(paste(myfolder,"/TMP_indlist.map",sep=""))
	Chip_ped <- read(paste(myfolder,"/TMP_indlist.ped",sep=""))
	rownames(Chip_ped) = Chip_ped[,2]; Chip_genotype = Chip_ped[,-c(1:6)]
	system(paste("rm -r ",myfolder,sep=""))

	ped_allele_merge <- function(k){
		x = Chip_ped[,(k*2+5):(k*2+6)]
		x_merge = paste(x[,1],x[,2],sep="")
		return(x_merge)
	}
	Chip_ped_merge = do.call('cbind',mclapply(1:((dim(Chip_ped)[2]-6)/2),ped_allele_merge,mc.cores=num_nodes))
	Chip_ped_merge = cbind(Chip_ped[,1:6],Chip_ped_merge)
	colnames(Chip_ped_merge) = c("FID","IID","PID","MID","sex","phe",as.character(Chip_map[,2]))
	
	allele_sum_all = c()
	for(j in 7:dim(Chip_ped_merge)[2]){
		allele_type = unique(Chip_ped_merge[,j])
		
		if(length(allele_type) == 3){
			allele_sum1 = length(Chip_ped_merge[,j][Chip_ped_merge[,j]==allele_type[1]])
			allele_sum2 = length(Chip_ped_merge[,j][Chip_ped_merge[,j]==allele_type[3]])
			allele_sum3 = length(Chip_ped_merge[,j][Chip_ped_merge[,j]==allele_type[3]])
			allele_sum = cbind(min(allele_sum1,allele_sum2,allele_sum3),colnames(Chip_ped_merge)[j])
			allele_sum_all = rbind(allele_sum_all,allele_sum)
		}
	}
	rownames(allele_sum_all) = allele_sum_all[,2]
	GT_sum_all[allele_sum_all[,2],i-3] = as.numeric(allele_sum_all[,1])
	
	print(i)
})
}

save(GT_sum_all,file="GT_Minimum_all.RData")

#######################################################################################
### STEP2. Draw the Minimum vs Sig-SNPs plot using Minimum from 15 to 50
#######################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)
setwd("/SAN/mottlab/heterosis/4_HSmice/3NonAdd_QTL/Venn-Plot")

load("/SAN/mottlab/heterosis/4_HSmice/3NonAdd_QTL/1_A+D/1_PheGen/0_Data/GT_Minimum_all.RData")
load("Venn_thrsuggest_data1.RData")
load("Venn_thrsuggest_data2.RData")

tmp_all = c()
for(i in 1:dim(GT_sum_all)[2]){
	tmp_name = paste(colnames(GT_sum_all)[i],rownames(GT_sum_all),sep=":")
	tmp = as.data.frame(GT_sum_all[,i]); rownames(tmp) = tmp_name
	tmp_all = rbind(tmp_all,tmp)
	print(i)
}

result_all_all = cbind(tmp_all[rownames(result_all_all),1],result_all_all)
save(result_all_all,file="Plot_GT_Minimum.RData")

result_all_all_clean = result_all_all[complete.cases(result_all_all),]; result_all_all_tmp = c()
for(i in seq(15,50,1)){
	result_all_all_tmp1 = subset(result_all_all_clean,result_all_all_clean[,1]>i-1 & result_all_all_clean[,"X.logP.DomCode."]>4.5)
	result_all_all_tmp2 = c(i,dim(result_all_all_tmp1)[1])
	result_all_all_tmp = rbind(result_all_all_tmp,result_all_all_tmp2)
}

png("Plot_GT_Minimum.png",width=10000,height=8000,res=600,type="cairo")
par(mai=c(1.5,1.6,0.7,0.7),mgp=c(4.7,2,0))
plot(result_all_all_tmp[,1],result_all_all_tmp[,2],col="red",xlab="Minimum of effective individuals with AA/AB/BB",ylab="Sum of significant loci across all traits",cex.axis=2.3,cex.lab=2.5,cex.main=3,cex=3,pch=16)
lines(result_all_all_tmp[,1],result_all_all_tmp[,2],col="black",type="l",lwd=4,lty=1)
dev.off()


