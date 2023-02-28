
### STEP1. Extract SNP names from bim file based on KEY_SNPs.txt (chr+pos) ###
rm(list=ls());options(stringsAsFactors=FALSE)

KEY_SNPs = read.table("KEY_SNPs.txt",header=T)
bed_file = read.table("HSmice.bim",header=F)

KEY_SNPs_name = c()
for(i in 1:nrow(KEY_SNPs)){
try({
	tmp_name = subset(bed_file,bed_file[,1]==KEY_SNPs[i,1] & bed_file[,4]==KEY_SNPs[i,2])[1,2]
	KEY_SNPs_name = c(KEY_SNPs_name, tmp_name)
	print(i)
})
}

write.table(as.data.frame(KEY_SNPs_name), file="KEY_SNPs_name.txt", row.names=F, col.names=F, quote=F)

### STEP2. Convert the tped format into anno format ###
plink --bfile HSmice --extract KEY_SNPs_name.txt --recode --transpose --out tmpfile
#cut -d ' ' -f 1-100 tmpfile.tped > tmpfile1.tped # In cast the tped is too big 
tped2vepAnno -i tmpfile.tped -o KEY_SNPs.anno
rm tmpfile*

### STEP3. Convert KEY_SNPs.txt (chr+pos) into LiftOver format (bed format)  ###
cat KEY_SNPs.txt | awk -v OFS='' '{print "chr",$1,":",$2,"-",$2}' | sed '1d' > KEY_SNPs_OLD.bed
# Use the online version of LiftOver to convert the SNP from an old ref-version into a new ref-version #


