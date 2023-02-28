
rm(list=ls());options(stringsAsFactors=FALSE)

datdir = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE")
group_name = "HSrats_amygdala"
cov_num = 3; threadnum = 10

outdir = paste0("/SAN/mottlab/heterosis/3_HSrats/5VC_gene/", group_name, "_Trpts/1data")
setwd(outdir)

### STEP1. Generate the mgrm information file ###
mgrm_info = rbind(paste0("VC_", group_name, "_add"), paste0("VC_", group_name, "_dom.d"))
write.table(mgrm_info, file=paste0(group_name,"_add_domi.txt"), row.names=F, col.names=F, quote=F)

### STEP2. Generate the genotype files and kinship realted files ###
system(paste0("plink --noweb --bfile ", datdir, "/", group_name, " --maf 0.05 --geno 0.1 --make-bed --out VC_", group_name))

system(paste0("gcta64 --bfile VC_", group_name, " --thread-num ", threadnum, " --autosome --make-grm-gz --out VC_", group_name, "_add"))
system(paste0("gcta64 --bfile VC_", group_name, " --thread-num ", threadnum, " --autosome --make-grm-d-gz --out VC_", group_name, "_dom"))

system(paste0("gcta64 --bfile VC_", group_name, " --thread-num ", threadnum, " --autosome --make-grm --out VC_", group_name, "_add"))
system(paste0("gcta64 --bfile VC_", group_name, " --thread-num ", threadnum, " --autosome --make-grm-d --out VC_", group_name, "_dom"))

### STEP3. Generate the phenotype file ###
rawphe = read.table(paste0(datdir, "/", group_name, ".phe"), header=T)
rawfam = read.table(paste0(datdir, "/", group_name, ".fam"), header=F)
outphe = cbind("famid"=rawfam[,1], rawphe[,-c(2:(cov_num+1))])
write.table(outphe, file=paste0("VC_", group_name, ".pheno"), row.names=F, col.names=F, quote=F)

### STEP4. Calculate the variance components ###
for(j in 1:(ncol(outphe)-2)){

    phe_name = colnames(outphe)[j+2]
    system(paste0("gcta64 --reml --mgrm ", group_name, "_add_domi.txt --pheno VC_", group_name, ".pheno --mpheno ", j, " --out ../2results/vc_", phe_name))
    print(j)

}

