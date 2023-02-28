
##############################################################################################################
### GWAS_A+D : ( Take Pig F2 data as an example )
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
source("/home/ucbtlcu/bin/scripts/1_GWAS_AddDom_QC.r"); 
source("/home/ucbtlcu/bin/scripts/2_GWAS_AddDom_Pvalue.r"); 
source("../3_GWAS_AddDom_Plot_All_AD.r"); # Plotting Traits >thrsuggest & based on AvsAD-Pvalues # 
#source("../3_GWAS_AddDom_Plot_All_D.r"); # Plotting Traits >thrgenome & based on D-Pvalues #
source("/home/ucbtlcu/bin/scripts/plotGWAS_pig.r"); 
source("/home/ucbtlcu/bin/scripts/plotRegion_pig.r")

# [2] Specify directory and covariates types variable and focused phenotypes list (if want run all phenotypes, just set "PheList_Choose=F") # 
indir = "/SAN/mottlab/heterosis/1_F2pigs/1data/RAW_DGE"
outdir = "/SAN/mottlab/heterosis/1_F2pigs/3NonAdd_eQTL/1_A+D_M"

covariates_types = c("f","f","n","n","n","n","n","n","n","n","n","n")
names(covariates_types) = c("sex","batch","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PheList = c("id","sex","batch","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",gsub("EMMA_Plot_","",gsub(".png","",dir("../Sig-M/"))))

# [3] Use three functions one by one #
#GWAS_AddDom_QC(indir=indir, outdir=outdir, Input_name="F2_DGE_M", Input_type="PLINK", Kinship_type="GCTA_ad", PheList_Choose=T, PheList=PheList, Phe_ResDone = F, Phe_NormDone = F, Normal_method = "QUANTILE", covariates_sum=12, Phe_IndMinimum = 200, Phe_Extreme = 5, GT_maf = 0.05, GT_missing = 0.1, num_nodes=15)

#GWAS_AddDom_Pvalue(indir=indir, outdir=outdir, Input_name="F2_DGE_M", Kinship_type="GCTA_ad", VarComponent_Method="GCTA_ad", PheList_Choose=F, PheList=PheList, Run_separated = F, Run_separated_length = 300000,covariates_sum=12, covariates_types=covariates_types, Phe_IndMinimum=200, GT_IndMinimum=10, num_nodes=15)

GWAS_AddDom_Plot(outdir=outdir, PheList_Choose=F, PheList=PheList, covariates_sum=12, RegionMan_chr_whole=F, RegionMan_chr_region = 15000000, num_nodes=15)

