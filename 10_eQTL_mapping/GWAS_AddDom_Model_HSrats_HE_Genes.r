
##############################################################################################################
### ADDO_AddDom :
### 1. PLINK Input Format: (1) file.tped & file.tfam (2) file.phe (1st column name should be "id"; From 2nd column should start from covariates 
### and the sex column should coded as 0=female and 1=male; The sex column is required!) (3) file.covs (1st column is phenotype name; 2nd 
### column is corresponding covariates and all covariates should be separated by ",")
### 2. GenABEL Input Format: file.ABEL.dat (Just contain one GenABEL type variable named "dat") & file.covs (1st column is phenotype name; 2nd 
### column is corresponding covariates and all covariates should be separated by ",")
### 3. Required Softwares: plink (v.1.90) & emmax-kin/gemma/gcta64 (just install needed one)
##############################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

# [1] Load library and function #
library(data.table); library(parallel); library(BEDMatrix)
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/ADDO_AddDom1_QC.r");
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/ADDO_AddDom2_Pvalue.r");
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/ADDO_AddDom3_Plot.r");
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/ADDO_AddDom4_IntePlot.r");
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/plotGWAS.r");
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/plotRegion.r")
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/plotGWAS_ddinfo.r")
source("/SAN/mottlab/heterosis/bin/Heterosis_Models/plotRegion_ddinfo.r")

# [2] Specify directory and covariates types variable and focused phenotypes list (if want run all phenotypes, just set "PheList_Choose=F") # 
indir = "/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes"
outdir = "/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL_Gene/1_A+D_heart"

covariates_types = c("n", "n", "n"); names(covariates_types) = c("age", "sex", "BW_at_day28_pi")
PheList = c("id", "age", "sex", "BW_at_day28_pi", "ENSRNOT00000000008", " ENSRNOT00000000009", "ENSRNOT00000000010")

# [3] Use three functions one by one #
#ADDO_AddDom1_QC(indir=indir, outdir=outdir, Input_name="HSrats_heart", Input_type="PLINK", Kinship_type="GCTA_ad", PheList_Choose=F, PheList=PheList, Phe_ResDone = F, Phe_NormDone = T, Normal_method = "QUANTILE", covariates_sum=3, covariates_types=covariates_types, Phe_IndMinimum = 100, Phe_Extreme = 5, GT_maf = 0.05, GT_missing = 0.1, num_nodes=15)

ADDO_AddDom2_Pvalue(indir=indir, outdir=outdir, Input_name="HSrats_heart", Kinship_type="GCTA_ad", VarComponent_Method="GCTA_ad", PheList_Choose=F, PheList=PheList, covariates_sum=3, Phe_IndMinimum=100, GT_IndMinimum=10, matrix_acceleration=F, num_nodes=4)

#ADDO_AddDom3_Plot(outdir=outdir, PheList_Choose=F, PheList=PheList, covariates_sum=3, RegionMan_chr_whole=F, RegionMan_chr_region = 2000000, chrs_sum=21, num_nodes=15)

#ADDO_AddDom4_IntePlot(outdir=outdir, PheList_Choose=F, PheList=PheList, covariates_sum=3, RegionMan_chr_whole=F, RegionMan_chr_region = 2000000, Plot_model="NvsAD", chrs_sum=21, num_nodes=15)

