
rm(list=ls());options(stringsAsFactors=FALSE)

library(data.table)
read <- function(...) as.data.frame(fread(header=F,colClasses="double",...))
read.header <- function(...) as.data.frame(fread(header=T,colClasses="double",...))

gene_1trpt = read.table("/SAN/mottlab/heterosis/3_HSrats/1data/KEY_Transcripts_1isoform.txt", header=T)

trpts_info = read.header("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE/HSrats_amygdala.phe")
# trpts_info = read.header("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE/HSrats_heart.phe")

gene_info = read.header("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes/HSrats_amygdala.phe")
# gene_info = read.header("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes/HSrats_heart.phe")

rownames(trpts_info) = paste0("rat", trpts_info[,1])
rownames(gene_info) = paste0("rat", gene_info[,1])

common_sample = intersect(rownames(trpts_info), rownames(gene_info))
trpts_common = trpts_info[common_sample, gene_1trpt[,1]]
gene_common = gene_info[common_sample, gene_1trpt[,2]]

cor_res = c()
for(i in 1:ncol(trpts_common)){
	tmp1 = cor.test(trpts_common[,i], gene_common[,i])$estimate
	tmp2 = cor.test(trpts_common[,i], gene_common[,i])$p.value
	cor_tmp = cbind(tmp1, tmp2)
	cor_res = rbind(cor_res, cor_tmp)
}
rownames(cor_res) = colnames(trpts_common)
colnames(cor_res) = c("cor", "pvalue")

save(trpts_common, gene_common, cor_res, file="Cor_gene_trpt_AM.Rdata")
# save(trpts_common, gene_common, cor_res, file="Cor_gene_trpt_HE.Rdata")
