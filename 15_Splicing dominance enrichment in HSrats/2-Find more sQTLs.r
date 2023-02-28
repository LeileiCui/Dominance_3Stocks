
rm(list=ls());options(stringsAsFactors=FALSE)

pos_info = read.table("/SAN/mottlab/heterosis/3_HSrats/1data/KEY_Transcripts_Genes.txt", header=T)

trpts_info = read.table("/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL/1_A+D_0.999999r2_Amygdala/SNP_Peak/PeakSNP_AD_thrsuggest.txt", header=T)
gene_info = read.table("/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL_Gene/1_A+D_amygdala/SNP_Peak/PeakSNP_AD_thrsuggest.txt", header=T)

# trpts_info = read.table("/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL/1_A+D_0.999999r2_Heart/SNP_Peak/PeakSNP_AD_thrsuggest.txt", header=T)
# gene_info = read.table("/SAN/mottlab/heterosis/3_HSrats/3NonAdd_eQTL_Gene/1_A+D_heart/SNP_Peak/PeakSNP_AD_thrsuggest.txt", header=T)


if(T){

	gene_1trpts = names(table(pos_info[, "gene_id"]))[table(pos_info[, "gene_id"])==1] # To get the gene list for genes has >1 transcripts #
	gene_1trpts_info = c()

	for(i in 1:length(gene_1trpts)){
		tmp_trpts = subset(pos_info, pos_info[, "gene_id"]==gene_1trpts[i])
		tmp_info = tmp_trpts[1,-2]
		gene_1trpts_info = rbind(gene_1trpts_info, tmp_info)
	}

	gene_1trpts_info = cbind("gene_id"=gene_1trpts, gene_1trpts_info)
	colnames(gene_1trpts_info) = c("G_id", "T1_id", "T1_ph_chr", "T1_ph_start", "T1_ph_end")
	rownames(gene_1trpts_info) = gene_1trpts_info[,1]

	gene_1trpts_res = cbind("G_id"=gene_1trpts_info[,1], "G_chr"=NA, "G_pos"=NA, "G_logP_AD"=NA, "G_OD"=NA, gene_1trpts_info[,2:5], "T1_chr"=NA, "T1_pos"=NA, "T1_logP_AD"=NA, "T1_OD"=NA)

	for(j in 1:nrow(gene_1trpts_res)){
		tmp_g = subset(gene_info, gene_info[,1]==gene_1trpts_res[j, "G_id"])
		tmp_t1 = subset(trpts_info, trpts_info[,1]==gene_1trpts_res[j, "T1_id"])
		tmp_t2 = subset(trpts_info, trpts_info[,1]==gene_1trpts_res[j, "T2_id"])

		if(nrow(tmp_g)>0) { gene_1trpts_res[gene_1trpts_res[j,"G_id"], c("G_chr", "G_pos", "G_logP_AD", "G_OD")] = tmp_g[1, c("chr", "pos", "logP_AD", "OD.A")]}
		if(nrow(tmp_t1)>0) { gene_1trpts_res[gene_1trpts_res[j,"G_id"], c("T1_chr", "T1_pos", "T1_logP_AD", "T1_OD")] = tmp_t1[1, c("chr", "pos", "logP_AD", "OD.A")]}
		if(nrow(tmp_t2)>0) { gene_1trpts_res[gene_1trpts_res[j,"G_id"], c("T2_chr", "T2_pos", "T2_logP_AD", "T2_OD")] = tmp_t2[1, c("chr", "pos", "logP_AD", "OD.A")]}

		print(j)
	}

	write.table(gene_1trpts_res, file="Amygdala_1trpts.txt", row.names=F, col.names=T, quote=F)

}


if(T){

	gene_2trpts = names(table(pos_info[, "gene_id"]))[table(pos_info[, "gene_id"])==2] # To get the gene list for genes has 2 transcripts #
	gene_2trpts_info = c()

	for(i in 1:length(gene_2trpts)){
		tmp_trpts = subset(pos_info, pos_info[, "gene_id"]==gene_2trpts[i])
		tmp_info = cbind(tmp_trpts[1,-2], tmp_trpts[2,-2])
		gene_2trpts_info = rbind(gene_2trpts_info, tmp_info)
	}

	gene_2trpts_info = cbind("gene_id"=gene_2trpts, gene_2trpts_info)
	colnames(gene_2trpts_info) = c("G_id", "T1_id", "T1_ph_chr", "T1_ph_start", "T1_ph_end", "T2_id", "T2_ph_chr", "T2_ph_start", "T2_ph_end")
	rownames(gene_2trpts_info) = gene_2trpts_info[,1]

	gene_2trpts_res = cbind("G_id"=gene_2trpts_info[,1], "G_chr"=NA, "G_pos"=NA, "G_logP_AD"=NA, "G_OD"=NA, gene_2trpts_info[,2:5], "T1_chr"=NA, "T1_pos"=NA, "T1_logP_AD"=NA, "T1_OD"=NA, gene_2trpts_info[,6:9], "T2_chr"=NA, "T2_pos"=NA, "T2_logP_AD"=NA, "T2_OD"=NA)

	for(j in 1:nrow(gene_2trpts_res)){
		tmp_g = subset(gene_info, gene_info[,1]==gene_2trpts_res[j, "G_id"])
		tmp_t1 = subset(trpts_info, trpts_info[,1]==gene_2trpts_res[j, "T1_id"])
		tmp_t2 = subset(trpts_info, trpts_info[,1]==gene_2trpts_res[j, "T2_id"])

		if(nrow(tmp_g)>0) { gene_2trpts_res[gene_2trpts_res[j,"G_id"], c("G_chr", "G_pos", "G_logP_AD", "G_OD")] = tmp_g[1, c("chr", "pos", "logP_AD", "OD.A")]}
		if(nrow(tmp_t1)>0) { gene_2trpts_res[gene_2trpts_res[j,"G_id"], c("T1_chr", "T1_pos", "T1_logP_AD", "T1_OD")] = tmp_t1[1, c("chr", "pos", "logP_AD", "OD.A")]}
		if(nrow(tmp_t2)>0) { gene_2trpts_res[gene_2trpts_res[j,"G_id"], c("T2_chr", "T2_pos", "T2_logP_AD", "T2_OD")] = tmp_t2[1, c("chr", "pos", "logP_AD", "OD.A")]}

		print(j)
	}

	write.table(gene_2trpts_res, file="Amygdala_2trpts.txt", row.names=F, col.names=T, quote=F)

}


if(T){

	gene_3trpts = names(table(pos_info[, "gene_id"]))[table(pos_info[, "gene_id"])==3] # To get the gene list for genes has 2 transcripts #
	gene_3trpts_info = c()

	for(i in 1:length(gene_3trpts)){
		tmp_trpts = subset(pos_info, pos_info[, "gene_id"]==gene_3trpts[i])
		tmp_info = cbind(tmp_trpts[1,-2], tmp_trpts[2,-2], tmp_trpts[3,-2])
		gene_3trpts_info = rbind(gene_3trpts_info, tmp_info)
	}

	gene_3trpts_info = cbind("gene_id"=gene_3trpts, gene_3trpts_info)
	colnames(gene_3trpts_info) = c("G_id", "T1_id", "T1_ph_chr", "T1_ph_start", "T1_ph_end", "T2_id", "T2_ph_chr", "T2_ph_start", "T2_ph_end", "T3_id", "T3_ph_chr", "T3_ph_start", "T3_ph_end")
	rownames(gene_3trpts_info) = gene_3trpts_info[,1]

	gene_3trpts_res = cbind("G_id"=gene_3trpts_info[,1], "G_chr"=NA, "G_pos"=NA, "G_logP_AD"=NA, "G_OD"=NA, gene_3trpts_info[,2:5], "T1_chr"=NA, "T1_pos"=NA, "T1_logP_AD"=NA, "T1_OD"=NA, gene_3trpts_info[,6:9], "T2_chr"=NA, "T2_pos"=NA, "T2_logP_AD"=NA, "T2_OD"=NA, gene_3trpts_info[,10:13], "T3_chr"=NA, "T3_pos"=NA, "T3_logP_AD"=NA, "T3_OD"=NA)

	for(j in 1:nrow(gene_3trpts_res)){
		tmp_g = subset(gene_info, gene_info[,1]==gene_3trpts_res[j, "G_id"])
		tmp_t1 = subset(trpts_info, trpts_info[,1]==gene_3trpts_res[j, "T1_id"])
		tmp_t2 = subset(trpts_info, trpts_info[,1]==gene_3trpts_res[j, "T2_id"])
		tmp_t3 = subset(trpts_info, trpts_info[,1]==gene_3trpts_res[j, "T3_id"])

		if(nrow(tmp_g)>0) { gene_3trpts_res[gene_3trpts_res[j,"G_id"], c("G_chr", "G_pos", "G_logP_AD", "G_OD")] = tmp_g[1, c("chr", "pos", "logP_AD", "OD.A")]}
		if(nrow(tmp_t1)>0) { gene_3trpts_res[gene_3trpts_res[j,"G_id"], c("T1_chr", "T1_pos", "T1_logP_AD", "T1_OD")] = tmp_t1[1, c("chr", "pos", "logP_AD", "OD.A")]}
		if(nrow(tmp_t2)>0) { gene_3trpts_res[gene_3trpts_res[j,"G_id"], c("T2_chr", "T2_pos", "T2_logP_AD", "T2_OD")] = tmp_t2[1, c("chr", "pos", "logP_AD", "OD.A")]}
		if(nrow(tmp_t3)>0) { gene_3trpts_res[gene_3trpts_res[j,"G_id"], c("T3_chr", "T3_pos", "T3_logP_AD", "T3_OD")] = tmp_t3[1, c("chr", "pos", "logP_AD", "OD.A")]}

		print(j)
	}

	write.table(gene_3trpts_res, file="Amygdala_3trpts.txt", row.names=F, col.names=T, quote=F)

}


################################################################################################################################################

QTL_TypeSum = function(x){
	tmp = abs(x)
	res1 = c(length(tmp[tmp<0.2]), length(tmp[tmp>0.2 & tmp<0.8]), length(tmp[tmp>0.8 & tmp<1.2]), length(tmp[tmp>1.2]))
	names(res1) = c("AD", "PD", "CD", "OD")

	res2 = c(length(tmp[tmp<0.8]), length(tmp[tmp>0.8]))
	names(res2) = c("Add", "Dom")

	return(list(res1, res2))
}

test_1trpts = subset(gene_1trpts_res, !is.na(gene_1trpts_res[,"T1_chr"]))[,6:13]
num_1trpts = QTL_TypeSum(test_1trpts[,8])[[2]]
num_1trpts_apco = QTL_TypeSum(test_1trpts[,8])[[1]]

test_2trpts1 = subset(gene_2trpts_res, !is.na(gene_2trpts_res[,"T1_chr"]))[,6:13]
test_2trpts2 = subset(gene_2trpts_res, !is.na(gene_2trpts_res[,"T2_chr"]))[,14:21]
colnames(test_2trpts2) = colnames(test_2trpts1)
test_2trpts = rbind(test_2trpts1, test_2trpts2)
num_2trpts = QTL_TypeSum(test_2trpts[,8])[[2]]
num_2trpts_apco = QTL_TypeSum(test_2trpts[,8])[[1]]

test_3trpts1 = subset(gene_3trpts_res, !is.na(gene_3trpts_res[,"T1_chr"]))[,6:13]
test_3trpts2 = subset(gene_3trpts_res, !is.na(gene_3trpts_res[,"T2_chr"]))[,14:21]
test_3trpts3 = subset(gene_3trpts_res, !is.na(gene_3trpts_res[,"T3_chr"]))[,22:29]
colnames(test_3trpts2) = colnames(test_3trpts1)
colnames(test_3trpts3) = colnames(test_3trpts1)
test_3trpts = rbind(test_3trpts1, test_3trpts2, test_3trpts3)
num_3trpts = QTL_TypeSum(test_3trpts[,8])[[2]]
num_3trpts_apco = QTL_TypeSum(test_3trpts[,8])[[1]]

rbind(num_1trpts, num_2trpts); chisq.test(rbind(num_1trpts, num_2trpts))
rbind(num_1trpts, num_3trpts); chisq.test(rbind(num_1trpts, num_3trpts))
rbind(num_1trpts, num_2trpts + num_3trpts); chisq.test(rbind(num_1trpts, num_2trpts + num_3trpts))

rbind(num_1trpts_apco, num_2trpts_apco); chisq.test(rbind(num_1trpts_apco, num_2trpts_apco))
rbind(num_1trpts_apco, num_3trpts_apco); chisq.test(rbind(num_1trpts_apco, num_3trpts_apco))
rbind(num_1trpts_apco, num_2trpts_apco + num_3trpts_apco); chisq.test(rbind(num_1trpts_apco, num_2trpts_apco + num_3trpts_apco))



