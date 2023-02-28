
rm(list=ls());options(stringsAsFactors=FALSE)

outdir = "/SAN/mottlab/heterosis/4_HSmice/5VC_gene/HSmice_liver"; setwd(outdir)
resname = "HSmice_liver"

hsq_list = list.files(path=paste0(outdir, "/2results"), pattern=".hsq")

tmp_result_all = c()
for(i in 1:length(hsq_list)){
	tmp_hsq = read.table(paste0("2results/",hsq_list[i]), header=T, fill=T)
	rownames(tmp_hsq) = tmp_hsq[,1]
	tmp_result = c(hsq_list[i],as.numeric(as.character(tmp_hsq["V(G1)/Vp","Variance"])),as.numeric(as.character(tmp_hsq["V(G1)/Vp","SE"])),as.numeric(as.character(tmp_hsq["V(G2)/Vp","Variance"])),as.numeric(as.character(tmp_hsq["V(G2)/Vp","SE"])))
	tmp_result_all = rbind(tmp_result_all,tmp_result)
	print(paste(i," ",hsq_list[i],sep=""))
}

tmp_result_all[,1] = gsub(".hsq","",tmp_result_all[,1])
colnames(tmp_result_all) = c("Trait","Add_Var","Add_SE","Dom_Var","Dom_SE")

write.table(tmp_result_all, file=paste0(resname, "_AddDom_Variance.txt"), row.names=F, col.names=T, quote=F)

