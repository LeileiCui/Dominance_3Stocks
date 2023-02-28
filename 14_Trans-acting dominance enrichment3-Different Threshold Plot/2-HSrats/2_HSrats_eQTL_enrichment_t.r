rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/6-4-Dom Enrichment-Different Threshold Plot/2-HSrats")

eQTL_filenames = c("HSrats_amygdala_eQTL_AD_t.txt", "HSrats_heart_eQTL_AD_t.txt")
res_filenames = c("HSrats_amygdala", "HSrats_heart")

for(name_sort in 1:length(eQTL_filenames)){

	eQTL_data = read.table(file=eQTL_filenames[name_sort], header=T, sep="\t")
	eQTL_data_rm = eQTL_data[complete.cases(eQTL_data),]; rownames(eQTL_data_rm) = 1:nrow(eQTL_data_rm)
	res_name = res_filenames[name_sort]

	### Calculate the cis/trans ratio within add/dom eQTLs ###

	eQTL_thres = seq(from=3.5, to=15, by=0.01)
	eQTL_results = as.data.frame(matrix(0, nrow=length(eQTL_thres), ncol=7))
	colnames(eQTL_results) = c("eQTL_thres", "Num_all", "Num_add_cis", "Num_add_trans", "Num_dom_cis", "Num_dom_trans", "Fisher_pvalue")

	for(j in 1:length(eQTL_thres)){

		try({

				eQTL_thres_tmp = eQTL_thres[j]
				eQTL_data_tmp = subset(eQTL_data_rm, eQTL_data_rm[,"logP_AD"]>eQTL_thres_tmp)
				eQTL_data_addtmp = subset(eQTL_data_rm, eQTL_data_rm[,"logP_AD"]>eQTL_thres_tmp & eQTL_data_rm[,"eQTL_Type"]%in%c("additive", "partial-dominant"))
				eQTL_data_domtmp = subset(eQTL_data_rm, eQTL_data_rm[,"logP_AD"]>eQTL_thres_tmp & eQTL_data_rm[,"eQTL_Type"]%in%c("dominant", "overdominant"))

				add_cis_tmp = nrow(subset(eQTL_data_addtmp, eQTL_data_addtmp[,"cis_trans"]=="cis"))
				add_trans_tmp = nrow(subset(eQTL_data_addtmp, eQTL_data_addtmp[,"cis_trans"]=="trans"))
				dom_cis_tmp = nrow(subset(eQTL_data_domtmp, eQTL_data_domtmp[,"cis_trans"]=="cis"))
				dom_trans_tmp = nrow(subset(eQTL_data_domtmp, eQTL_data_domtmp[,"cis_trans"]=="trans"))

				pvalue_tmp = fisher.test(matrix(c(add_cis_tmp, dom_cis_tmp, add_trans_tmp, dom_trans_tmp), 2, 2))$p.value

			})

		if(nrow(eQTL_data_tmp)>0){
			
			eQTL_results[j,] = c(eQTL_thres_tmp, nrow(eQTL_data_tmp), add_cis_tmp, add_trans_tmp, dom_cis_tmp, dom_trans_tmp, pvalue_tmp)
		
		}

	}

	write.table(eQTL_results, file=paste0("Summary_", res_name, "_eQTL_enrichment_t.txt"), row.names=F, col.names=T, quote=F, sep="\t")

}

### Draw the dual-Y-axis plot for enrichment pvalue and eQTL sum ###

rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/6-4-Dom Enrichment-Different Threshold Plot/2-HSrats")

res_filenames = c("F2pigs_liver", "F2pigs_muscle", "HSrats_amygdala", "HSrats_heart", "HSmice_hippocampus", "HSmice_liver", "HSmice_lung")[3:4]
plot_cols = c("black", "bisque4", "bisque")[1:2]

ylim_max1_all = c(); ylim_max2_all = c()
for(name_sort in 1:length(res_filenames)){

	eQTL_results_tmp = read.table(paste0("Summary_", res_filenames[name_sort], "_eQTL_enrichment_t.txt"), header=T)

	data_x = eQTL_results_tmp[,"eQTL_thres"]
	data_y1 = -log10(eQTL_results_tmp[,"Fisher_pvalue"])
	data_y2 = eQTL_results_tmp[,"Num_all"]

	ylim_max1_tmp = max(data_y1[is.finite(data_y1)])*3; ylim_max2_tmp = max(data_y2[is.finite(data_y2)])*1.12
	ylim_max1_all = c(ylim_max1_all, ylim_max1_tmp); ylim_max2_all = c(ylim_max2_all, ylim_max2_tmp); 

}

for(name_sort in 1:length(res_filenames)){

	eQTL_results_tmp = read.table(paste0("Summary_", res_filenames[name_sort], "_eQTL_enrichment_t.txt"), header=T)

	data_x = eQTL_results_tmp[,"eQTL_thres"]
	data_y1 = -log10(eQTL_results_tmp[,"Fisher_pvalue"])
	data_y2 = eQTL_results_tmp[,"Num_all"]

	lwd_used=3; xlim_max = 15; ylim_max1 = max(ylim_max1_all); ylim_max2 = max(ylim_max2_all)

	if(name_sort == 1){

		png(paste0("Plot_eQTL_enrichment_t.png"), width=4300, height=4000, res=600, type="cairo")

		par(mar=c(4,4,2,4.5), mgp=c(2.2,0.8,0), cex.main=1.5)
		plot(data_x, data_y1, type="l", col=plot_cols[name_sort], xlab="Significant thresholds of eQTLs (HSrats_t)", ylab="-logP of Fisher test of eQTL trans-acting enrichment", lty=1, xlim=c(0,xlim_max), ylim=c(0,ylim_max1), lwd=lwd_used)

		par(mar=c(4,4,2,4.5), mgp=c(2.2,0.8,0), cex.main=1.5, new=T)
		plot(data_x, data_y2, type="l", col=plot_cols[name_sort], axes=F, xlab="", ylab="", lty=2, xlim=c(0,xlim_max), ylim=c(0,ylim_max2), lwd=lwd_used)
		axis(side = 4, at = pretty(c(0,ylim_max2)))

	}else{

		par(mar=c(4,4,2,4.5), mgp=c(2.2,0.8,0), cex.main=1.5, new=T)
		plot(data_x, data_y1, type="l", col=plot_cols[name_sort], xlab="", ylab="", lty=1, xlim=c(0,xlim_max), ylim=c(0,ylim_max1), lwd=lwd_used)

		par(mar=c(4,4,2,4.5), mgp=c(2.2,0.8,0), cex.main=1.5, new=T)
		plot(data_x, data_y2, type="l", col=plot_cols[name_sort], axes=F, xlab="", ylab="", lty=2, xlim=c(0,xlim_max), ylim=c(0,ylim_max2), lwd=lwd_used)
		axis(side = 4, at = pretty(c(0,ylim_max2)))

	}

	if(name_sort == length(res_filenames)){

		legend("topleft", cex=1, legend=res_filenames, col=plot_cols, lty=1, lwd=lwd_used)
		legend("topright", cex=1, legend=c("Left y-axis (Fisher_pvalue)","Right y-axis (Num_eQTLs)"), lty=c(1,2), lwd=lwd_used)

		mtext("Num of all significant eQTLs", side = 4, line = 2)
		dev.off()

	}

}
