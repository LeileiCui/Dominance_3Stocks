
##############################################################################################################
### GWAS_2_T-statistics : ( Take Rat HS data as an example )
### 1. PLINK Input Format: (1) file.tped & file.tfam (2) file.phe (1st column name should be "id"; From 2nd column should start from covariates 
### and the sex column should coded as 0=female and 1=male; The sex column is necessary!) (3) file.covs (1st column is phenotype name; 2nd 
### column is corresponding covariates and all covariates should be separated by ",")
### 2. GenABEL Input Format: file.ABEL.dat (Just contain one GenABEL type variable named "dat") & file.covs (1st column is phenotype name; 2nd 
### column is corresponding covariates and all covariates should be separated by ",")
### 3. Required Softwares: plink2 (v.1.90) & emmax-kin/gemma/gcta64 (just install needed one)
##############################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

# [1] Load library and function #
library(data.table); library(parallel); library(bigmemory); library(mvtnorm); library(MASS)
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/1_GWAS_Tstatistics_QC.r")
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/2_GWAS_Tstatistics_Pvalue.r")
source("/SAN/mottlab/heterosis/3_HSrats/3NonAdd_QTL/RawPhe_Plot_Script/3_GWAS_Tstatistics_Plot.r")
source("/SAN/mottlab/heterosis/bin/scripts/plotGWAS_rat.r")
source("/SAN/mottlab/heterosis/bin/scripts/plotRegion_rat.r")

# [2] Specify directory and covariates types variable and focused phenotypes list (if want run all phenotypes, just set "PheList_Choose=F") # 
indir = "/SAN/mottlab/heterosis/3_HSrats/1data"
outdir = "/SAN/mottlab/heterosis/3_HSrats/3NonAdd_QTL/2_T-statistics"

covariates_types = c("f","f","f","n","n","sel","n","sel","n","n")
names(covariates_types) = c("sex","batch","is.albino","Haemalysis","BW_at_day28_pi","test_worked","age","reliable","BW_at_IPGTT","BW_at_day9_pi")
PheList = c("id","age","batch","BW_at_day28_pi","BW_at_day9_pi","BW_at_IPGTT","Haemalysis","is.albino","reliable","sex","test_worked","ALP_bc","dist_femur_CRT_A","fem_neck_TOT_DEN","fem_neck_TRAB_A","lumbar_Area","NW","BD","died","CD4_CD8_ratio","CD8Tcells_Abs","pctCD25posCD4","Has1kidney","HDW_bc","MPM_bc","AA_IL_nb_pos_qn","AA_IL_nb_qn","AA_IL_score_pos_qn","AA_IL_score_qn","AA_nb_qn","AA_score_qn","IL_nb_qn","IL_score_qn","old_AA_IL_score_qn","Avoidances1_5_qn")

# [3] Use three functions one by one #
#GWAS_Tstatistics_QC(indir=indir, outdir=outdir, Input_name="HSrats", Input_type="PLINK", Kinship_type="GCTA_ad", PheList_Choose=F, PheList=PheList, Phe_ResDone = F, Phe_NormDone = T, Normal_method = "QUANTILE", covariates_sum=10, Phe_IndMinimum = 200, Phe_Extreme = 5, GT_maf = 0.05, GT_missing = 0.1, num_nodes=15)

GWAS_Tstatistics_Pvalue(indir=indir, outdir=outdir, Input_name="HSrats", Kinship_type="GCTA_ad", VarComponent_Method="GCTA_ad", PheList_Choose=F, PheList=PheList, Run_separated = F,covariates_sum=10, covariates_types=covariates_types, Phe_IndMinimum=200, GT_IndMinimum=10, num_nodes=30)

#GWAS_Tstatistics_Plot(outdir=outdir, PheList_Choose=F, PheList=PheList, covariates_sum=10, RegionMan_chr_whole=F, RegionMan_chr_region = 15000000, num_nodes=15)


