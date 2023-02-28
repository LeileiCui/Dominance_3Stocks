
##############################################################################################################
### GWAS_Tstatistics : ( Take Pig F2 data as an example )
### 1. PLINK Input Format: (1) file.tped & file.tfam (2) file.phe (1st column name should be "id"; From 2nd column should start from covariates 
### and the sex column should coded as 0=female and 1=male; The sex column is necessary!) (3) file.covs (1st column is phenotype name; 2nd 
### column is corresponding covariates and all covariates should be separated by ",")
### 2. GenABEL Input Format: file.ABEL.dat (Just contain one GenABEL type variable named "dat") & file.covs (1st column is phenotype name; 2nd 
### column is corresponding covariates and all covariates should be separated by ",")
### 3. Required Softwares: plink2 (v.1.90) & emmax-kin/gemma/gcta64 (just install needed one)
##############################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

# [1] Load library and function #
library(GenABEL); library(emma); library(data.table); library(parallel); library(mvtnorm)
source("/home/ucbtlcu/bin/scripts/1_GWAS_Tstatistics_QC.r"); 
source("/home/ucbtlcu/bin/scripts/2_GWAS_Tstatistics_Pvalue.r"); 
source("/home/ucbtlcu/bin/scripts/3_GWAS_Tstatistics_Plot.r"); 
source("/home/ucbtlcu/bin/scripts/plotGWAS_pig.r"); 
source("/home/ucbtlcu/bin/scripts/plotRegion_pig.r")

# [2] Specify directory and covariates types variable and focused phenotypes list (if want run all phenotypes, just set "PheList_Choose=F") # 
indir = "/SAN/mottlab/heterosis/1_F2pigs/1data"
outdir = "/SAN/mottlab/heterosis/1_F2pigs/3NonAdd_QTL/2_T-statistics"

covariates_types = c("f","f","n","n"); names(covariates_types) = c("sex","batch","age","bodyweight")
PheList = c("id","sex","batch","A14_0","A16_0","A18_3","A24_0","ADG120_240","ADG21_46","ADG46","An_3","AntebrachialL","aveBF","Birthweight","CruralL","D14_0","D16_0","D20_0","D20_1","D20_2","D20_5","Day46weight","fatcell_Area","fatcell_Perimeter","fatcell_V","foodfrequency","foodintake","GSP","HGB_18","LDa_45min","LDb_45min","LD_bag24drip","LDc_45min","LD_GLY_energy","LDH_45min","LD_lac","LDpH9_24","LDT3h","LDT45min","LeftEarArea","LeftEarW","LFLS213","LiverW","MCH_240","MCH_46","MCHC_46","MCV_18","MPV_240","NeckW","PLT_240","RBC_18","RBC_46","RDW_120","RDW_46","RFLS213","RightEarArea","RightEarW","SMa_24h","SM_bag24drip","SMc_45min","SMcolor_45min","SMH_24h","SML_24h","SMMYO","SMpH15h","SMpH24h","SMpH9_24","SMpH9h","SMT24h","SMT45min","SpleenW","WaistBF","WaistVN","WBC_120")

# [3] Use three functions one by one #
GWAS_Tstatistics_QC(indir=indir, outdir=outdir, Input_name="F2", Input_type="GenABEL", Kinship_type="GCTA_ad", PheList_Choose=F, PheList=PheList, Phe_ResDone = F, Phe_NormDone = F, Normal_method = "QUANTILE", covariates_sum=2, Phe_IndMinimum = 200, Phe_Extreme = 5, GT_maf = 0.05, GT_missing = 0.1, num_nodes=15)

GWAS_Tstatistics_Pvalue(indir=indir, outdir=outdir, Input_name="F2", Kinship_type="GCTA_ad", VarComponent_Method="GCTA_ad", PheList_Choose=F, PheList=PheList, Run_separated = F, Run_separated_length = 300000, covariates_sum=2, covariates_types=covariates_types, Phe_IndMinimum=200, GT_IndMinimum=10, num_nodes=15)

GWAS_Tstatistics_Plot(outdir=outdir, PheList_Choose=F, PheList=PheList, covariates_sum=2, RegionMan_chr_whole=F, RegionMan_chr_region = 15000000, num_nodes=15)

