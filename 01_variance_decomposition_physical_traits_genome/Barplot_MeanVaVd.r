

rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/2-全基因组方差剖分-3柱状图")


### [1] Seperated barplots of average Va and Vd of physical traits ###

qtl1 = read.table("QTL-1.txt", header=T, sep="\t")
qtl2 = read.table("QTL-2.txt", header=T, sep="\t")
qtl3 = read.table("QTL-3.txt", header=T, sep="\t")

qtl_plotdat = t(rbind( c(mean(qtl1[,3]), mean(qtl1[,5])), c(mean(qtl1[,7],na.rm=T), mean(qtl1[,8],na.rm=T)),
                       c(mean(qtl2[,3]), mean(qtl2[,5])), c(mean(qtl2[,7],na.rm=T), mean(qtl2[,8],na.rm=T)),
     				   c(mean(qtl3[,3]), mean(qtl3[,5])), c(mean(qtl3[,7],na.rm=T), mean(qtl3[,8],na.rm=T)) ))

rownames(qtl_plotdat) = c("Va", "Vd")
colnames(qtl_plotdat) = usedname =  c("F2 Pigs Whole-Genome", "F2 Pigs QTLs", "HS Rats Whole-Genome", "HS Rats QTLs", "HS Mice Whole-Genome", "HS Mice QTLs")

for(i in 1:6){

	tmp_qtl_plotdat = as.matrix(qtl_plotdat[,i]); colnames(tmp_qtl_plotdat) = usedname[i]

	png(paste0("1-", i, "-Barplot_MeanVaVd_PhyTraits.png"), width=2800, height=3300, res=600)

	par(mar=c(2,3,1.3,1), mgp=c(1.7,0.5,0), cex.main=0.9, lwd=1)
	barplot(tmp_qtl_plotdat, beside=T, col=c("#0000FFB3","#FF0000B3"), border="black", ylim=c(-0.07, 0.55), ylab="Phenotypic variance explained by additive (Va) and dominance (Vd)", cex.axis=1.2, cex.lab=1, cex.names=2.2, font=2)

	par(mar=c(2,3,1.3,1), mgp=c(1.7,0.5,0), cex.main=0.9, lwd=4.5)
	barplot(tmp_qtl_plotdat, beside=T, col=c("white","white"), border="black", ylim=c(-0.07, 0.55), ylab="Phenotypic variance explained by additive (Va) and dominance (Vd)", cex.axis=1.2, cex.lab=1, cex.names=2.2, angle = c(45,45), density = c(4,4), font=2, add=T)	

	text(paste0(round(as.numeric(qtl_plotdat[,i])*100, 2), "%"), x=c(1.5,2.5), y=as.numeric(qtl_plotdat[,i])+0.04, cex=2.7, col=c("#0000FFB3","#FF0000B3"), font=2)
	text(c("Va","Vd"), x=c(1.5,2.5), y=c(-0.03, -0.03), cex=1.8, col=c("#0000FFB3","#FF0000B3"), font=2)

	dev.off()

}


### [2] Integrated barplot of average Va and Vd of gene expression traits ###

eqtl1 = read.table("eQTL-1.txt", header=T)
eqtl2 = read.table("eQTL-2.txt", header=T)
eqtl3 = read.table("eQTL-3.txt", header=T)
eqtl4 = read.table("eQTL-4.txt", header=T)
eqtl5 = read.table("eQTL-5.txt", header=T)
eqtl6 = read.table("eQTL-6.txt", header=T)
eqtl7 = read.table("eQTL-7.txt", header=T)

eqtl_plotdat = t(rbind( c(mean(eqtl1[,2]), mean(eqtl1[,4])),
						c(mean(eqtl2[,2]), mean(eqtl2[,4])),
						c(mean(eqtl3[,2]), mean(eqtl3[,4])),
						c(mean(eqtl4[,2]), mean(eqtl4[,4])),
						c(mean(eqtl5[,2]), mean(eqtl5[,4])),
						c(mean(eqtl6[,2]), mean(eqtl6[,4])),
						c(mean(eqtl7[,2]), mean(eqtl7[,4])) ))

rownames(eqtl_plotdat) = c("Va", "Vd")
colnames(eqtl_plotdat) = c("F2 Pigs Liver","F2 Pigs Muscle","HS Rats Amygdala (Genes)","HS Rats Heart (Genes)","HS Mice Hippocampus","HS Mice Liver","HS Mice Lung")

png("2-Barplot_MeanVaVd_GeneTraits.png", width=2800*7, height=4400, res=600)

par(mar=c(2,3,1.3,1), mgp=c(1.7,0.5,0), cex.main=0.9, lwd=1)
barplot(as.matrix(eqtl_plotdat), beside=T, col=rep(c("#0000FFB3","#FF0000B3"),7), border="black", ylim=c(-0.014, 0.2), ylab="Phenotypic variance explained by additive (Va) and dominance (Vd)", cex.axis=1, cex.lab=1.4, cex.names=1.85, font=2)

par(mar=c(2,3,1.3,1), mgp=c(1.7,0.5,0), cex.main=0.9, lwd=4.5)
barplot(as.matrix(eqtl_plotdat), beside=T, col=rep(c("white","white"),7), border="black", ylim=c(-0.014, 0.2), ylab="Phenotypic variance explained by additive (Va) and dominance (Vd)", cex.axis=1, cex.lab=1.4, cex.names=1.85, angle = c(45,45), density = c(4,4), font=2, add=T)

tmpxpos = c(1.5,2.5, 4.5,5.5, 7.5,8.5, 10.5,11.5, 13.5,14.5, 16.5,17.5, 19.5,20.5)
text(paste0(round(as.numeric(eqtl_plotdat)*100, 2), "%"), x=tmpxpos, y=as.numeric(eqtl_plotdat)+0.01, cex=2.5, col=c("#0000FFB3","#FF0000B3"), font=2)
text(rep(c("Va","Vd"),7), x=tmpxpos, y=rep(-0.006, 7*2), cex=1.7, col=c("#0000FFB3","#FF0000B3"), font=2)

dev.off()


### [3] Seperated barplots of average Va and Vd of physical traits ###

qtl1 = read.table("QTL-1.txt", header=T, sep="\t")
qtl2 = read.table("QTL-2.txt", header=T, sep="\t")
qtl3 = read.table("QTL-3.txt", header=T, sep="\t")

col1 = read.table("1_F2pigs_col.txt", header=T, sep="\t")
col2 = read.table("2_HSrats_col.txt", header=T, sep="\t")
col3 = read.table("3_HSmice_col.txt", header=T, sep="\t")

col_dic = read.table("0_col_dic.txt", header=F, sep="\t", comment.char = "$"); rownames(col_dic) = col_dic[,1]

# [3-1] Pig Plot #

uniq_traitclass1 = unique(qtl1[,2]); res_traitclass1 = c()
for(i in 1:length(uniq_traitclass1)){
	tmp_dat = subset(qtl1, qtl1[,2]==uniq_traitclass1[i])
	tmp_res = c(mean(tmp_dat[,3]), mean(tmp_dat[,5]), mean(tmp_dat[,7],na.rm=T), mean(tmp_dat[,8],na.rm=T))
	res_traitclass1 = rbind(res_traitclass1, tmp_res)
}
rownames(res_traitclass1) = uniq_traitclass1; colnames(res_traitclass1) = c("Va", "Vd", "Va_qtl", "Vd_qtl")
# res_traitclass1 = res_traitclass1[order(res_traitclass1[,1],decreasing=T),] # order by Va #
res_traitclass1 = res_traitclass1[order(res_traitclass1[,2],decreasing=T),] # order by Vd #

uniq_col1 = unique(col1[,2:3]); rownames(uniq_col1) = uniq_col1[,1]
res_traitclass1 = cbind(res_traitclass1, "col" = uniq_col1[rownames(res_traitclass1),2])
plotdat1 = cbind(res_traitclass1[,1:2], col_dic[res_traitclass1[,"col"],c(5,6)])

for(i in 1:9){

	plotdat1[,1] = as.numeric(plotdat1[,1]); plotdat1[,2] = as.numeric(plotdat1[,2]); plotdat = plotdat1

	if(i == 1){ 
		
		png("3-1-Barplot_MeanVaVd_PhyTraits.png", width=2800*3, height=3300*3, res=600)
		layout(t(matrix(1:9, nrow=3)), widths=rep(1,3), heights=rep(1,3))
	}

	par(mar=c(2,3,1.3,1), mgp=c(1.7,0.5,0), cex.main=0.9)
	tmp_qtl_plotdat = t(as.matrix(plotdat[i,1:2]))
	barplot(tmp_qtl_plotdat, beside=T, col=as.character(plotdat[i,3:4]), border="black", ylim=c(-0.07, 0.55), ylab="Phenotypic variance explained by additive (Va) and dominance (Vd)", cex.axis=1.2, cex.lab=1, cex.names=3.4, font=2)

	text(paste0(round(as.numeric(plotdat[i,1:2])*100, 2), "%"), x=c(1.5,2.5), y=as.numeric(plotdat[i,1:2])+0.04, cex=3.8, col=c("blue","red"), font=2)
	text(c("Va","Vd"), x=c(1.5,2.5), y=c(-0.03, -0.03), cex=2.4, col=c("blue","red"), font=2)

	if(i == 9){ dev.off() }

}


# [3-2] Rat Plot #

uniq_traitclass2 = unique(qtl2[,2]); res_traitclass2 = c()
for(i in 1:length(uniq_traitclass2)){
	tmp_dat = subset(qtl2, qtl2[,2]==uniq_traitclass2[i])
	tmp_res = c(mean(tmp_dat[,3]), mean(tmp_dat[,5]), mean(tmp_dat[,7],na.rm=T), mean(tmp_dat[,8],na.rm=T))
	res_traitclass2 = rbind(res_traitclass2, tmp_res)
}
rownames(res_traitclass2) = uniq_traitclass2; colnames(res_traitclass2) = c("Va", "Vd", "Va_qtl", "Vd_qtl")
# res_traitclass2 = res_traitclass2[order(res_traitclass2[,1],decreasing=T),] # order by Va #
res_traitclass2 = res_traitclass2[order(res_traitclass2[,2],decreasing=T),] # order by Vd #

uniq_col2 = unique(col2[,2:3]); rownames(uniq_col2) = uniq_col2[,1]
res_traitclass2 = cbind(res_traitclass2, "col" = uniq_col2[rownames(res_traitclass2),2])
plotdat2 = cbind(res_traitclass2[,1:2], col_dic[res_traitclass2[,"col"],c(5,6)])

for(i in 1:12){

	plotdat2[,1] = as.numeric(plotdat2[,1]); plotdat2[,2] = as.numeric(plotdat2[,2]); plotdat = plotdat2[-13,]
	rownames(plotdat)[5] = "Cardiovascular function"; rownames(plotdat)[8] = "Induced neuroinflammation"; rownames(plotdat)[12] = "Bone mass and strength"

	if(i == 1){ 
		
		png("3-2-Barplot_MeanVaVd_PhyTraits.png", width=2800*3, height=3300*3, res=600)
		layout(t(matrix(1:12, nrow=4)), widths=rep(1,4), heights=rep(1,3))
	}

	par(mar=c(2,3,1.5,1), mgp=c(1.7,0.5,0), cex.main=0.9)
	tmp_qtl_plotdat = t(as.matrix(plotdat[i,1:2]))
	barplot(tmp_qtl_plotdat, beside=T, col=as.character(plotdat[i,3:4]), border="black", ylim=c(-0.07, 0.69), ylab="Phenotypic variance explained by additive (Va) and dominance (Vd)", cex.axis=1.2, cex.lab=1, cex.names=2.4, font=2)

	text(paste0(round(as.numeric(plotdat[i,1:2])*100, 2), "%"), x=c(1.5,2.5), y=as.numeric(plotdat[i,1:2])+0.04, cex=3.5, col=c("blue","red"), font=2)
	text(c("Va","Vd"), x=c(1.5,2.5), y=c(-0.03, -0.03), cex=2.4, col=c("blue","red"), font=2)

	if(i == 12){ dev.off() }

}


# [3-3] Mouse Plot #

uniq_traitclass3 = unique(qtl3[,2]); res_traitclass3 = c()
for(i in 1:length(uniq_traitclass3)){
	tmp_dat = subset(qtl3, qtl3[,2]==uniq_traitclass3[i])
	tmp_res = c(mean(tmp_dat[,3]), mean(tmp_dat[,5]), mean(tmp_dat[,7],na.rm=T), mean(tmp_dat[,8],na.rm=T))
	res_traitclass3 = rbind(res_traitclass3, tmp_res)
}
rownames(res_traitclass3) = uniq_traitclass3; colnames(res_traitclass3) = c("Va", "Vd", "Va_qtl", "Vd_qtl")
# res_traitclass3 = res_traitclass3[order(res_traitclass3[,1],decreasing=T),] # order by Va #
res_traitclass3 = res_traitclass3[order(res_traitclass3[,2],decreasing=T),] # order by Vd #

uniq_col3 = unique(col3[,2:3]); rownames(uniq_col3) = uniq_col3[,1]
res_traitclass3 = cbind(res_traitclass3, "col" = uniq_col3[rownames(res_traitclass3),2])
plotdat3 = cbind(res_traitclass3[,1:2], col_dic[res_traitclass3[,"col"],c(5,6)])


for(i in 1:12){

	plotdat3[,1] = as.numeric(plotdat3[,1]); plotdat3[,2] = as.numeric(plotdat3[,2]); plotdat = plotdat3[-c(6,14,15),]

	if(i == 1){ 
		
		png("3-3-Barplot_MeanVaVd_PhyTraits.png", width=2800*3, height=3300*3, res=600)
		layout(t(matrix(1:12, nrow=4)), widths=rep(1,4), heights=rep(1,3))
	}

	par(mar=c(2,3,1.5,1), mgp=c(1.7,0.5,0), cex.main=0.9)
	tmp_qtl_plotdat = t(as.matrix(plotdat[i,1:2]))
	barplot(tmp_qtl_plotdat, beside=T, col=as.character(plotdat[i,3:4]), border="black", ylim=c(-0.07, 0.45), ylab="Phenotypic variance explained by additive (Va) and dominance (Vd)", cex.axis=1.2, cex.lab=1, cex.names=2.4, font=2)

	text(paste0(round(as.numeric(plotdat[i,1:2])*100, 2), "%"), x=c(1.5,2.5), y=as.numeric(plotdat[i,1:2])+0.04, cex=3.5, col=c("blue","red"), font=2)
	text(c("Va","Vd"), x=c(1.5,2.5), y=c(-0.03, -0.03), cex=2.4, col=c("blue","red"), font=2)

	if(i == 12){ dev.off() }

}






