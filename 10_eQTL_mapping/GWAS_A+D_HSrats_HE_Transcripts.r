
##############################################################################################################
### GWAS_A+D : ( Take Rat HS data as an example )
### 1. PLINK Input Format: (1) file.tped & file.tfam (2) file.phe (1st column name should be "id"; From 2nd column should start from covariates 
### and the sex column should coded as 0=female and 1=male; The sex column is necessary!) (3) file.covs (1st column is phenotype name; 2nd 
### column is corresponding covariates and all covariates should be separated by ",")
### 2. GenABEL Input Format: file.ABEL.dat (Just contain one GenABEL type variable named "dat") & file.covs (1st column is phenotype name; 2nd 
### column is corresponding covariates and all covariates should be separated by ",")
### 3. Required Softwares: plink2 (v.1.90) & emmax-kin/gemma/gcta64 (just install needed one)
##############################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

# [1] Load library and function #
library(GenABEL); library(emma); library(data.table); library(parallel)
source("/SAN/mottlab/heterosis/bin/scripts/1_GWAS_AddDom_QC.r"); 
source("/SAN/mottlab/heterosis/bin/scripts/2_GWAS_AddDom_Pvalue.r"); 
source("../3_GWAS_AddDom_Plot_All_AD.r"); # Plotting Traits >thrsuggest & based on NvsAD-Pvalues (All in 1 Plot) # 
#source("../3_GWAS_AddDom_Plot_Sig_AD.r"); # Plotting Traits >thrsuggest & based on NvsAD-Pvalues (Separated Plots) # 
#source("../3_GWAS_AddDom_Plot_All_D.r"); # Plotting Traits >thrsuggest & based on D-Pvalues (All in 1 Plot) # 
#source("../3_GWAS_AddDom_Plot_Sig_D.r"); # Plotting Traits >thrsuggest & based on D-Pvalues (Separated Plots) # 
source("/SAN/mottlab/heterosis/bin/scripts/plotGWAS_rat.r"); 
source("/SAN/mottlab/heterosis/bin/scripts/plotRegion_rat.r")

# [2] Specify directory and covariates types variable and focused phenotypes list (if want run all phenotypes, just set "PheList_Choose=F") # 
indir = "/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE"
outdir = "/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL/1_A+D_0.999999r2_Heart"

covariates_types = c("n","f","n")
names(covariates_types) = c("age","sex","BW_at_day28_pi")
PheList = c("id","age","sex","BW_at_day28_pi","ENSRNOT00000037224","ENSRNOT00000008936","ENSRNOT00000006604","ENSRNOT00000061800","ENSRNOT00000014171","ENSRNOT00000061799","ENSRNOT00000006157","ENSRNOT00000006151","ENSRNOT00000011307","ENSRNOT00000016193","ENSRNOT00000022823","ENSRNOT00000000478","ENSRNOT00000000497","ENSRNOT00000000641","ENSRNOT00000065629")

# [3] Use three functions one by one #
#GWAS_AddDom_QC(indir=indir, outdir=outdir, Input_name="HSrats_0.999999r2_heart", Input_type="PLINK", Kinship_type="GCTA_ad", PheList_Choose=F, PheList=PheList, Phe_ResDone = F, Phe_NormDone = T, Normal_method = "QUANTILE", covariates_sum=3, Phe_IndMinimum = 100, Phe_Extreme = 5, GT_maf = 0.05, GT_missing = 0.1, num_nodes=15)

#GWAS_AddDom_Pvalue(indir=indir, outdir=outdir, Input_name="HSrats_0.999999r2_heart", Kinship_type="GCTA_ad", VarComponent_Method="GCTA_ad", PheList_Choose=F, PheList=PheList, Run_separated = F, Run_separated_length = 300000,covariates_sum=3, covariates_types=covariates_types, Phe_IndMinimum=100, GT_IndMinimum=10, num_nodes=15)

GWAS_AddDom_Plot(outdir=outdir, PheList_Choose=F, PheList=PheList, covariates_sum=3, RegionMan_chr_whole=F, RegionMan_chr_region = 15000000, num_nodes=15)

