
######################################################################################################################
### Plot1: Scatter plot of all individuals-pairs using off-diagonal elements of additive GRM & dominance GRM
######################################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

outdir = "/SAN/mottlab/heterosis/3_HSrats/5VC_gene/HSrats_heart_Trpts"; setwd(outdir)
resname = "HSrats_heart"

grmAdd = read.table(file=paste0("1data/VC_", resname, "_add.grm.gz"))
grmDom = read.table(file=paste0("1data/VC_", resname, "_dom.d.grm.gz"))
grmAdd = grmAdd[grmAdd[,1]!=grmAdd[,2],]
grmDom = grmDom[grmDom[,1]!=grmDom[,2],]

res_grm = cbind(indpair=paste(grmAdd[,1],"_",grmAdd[,2],sep=""),grmAdd=grmAdd[,4],grmDom=grmDom[,4])
write.table(res_grm,file=paste0(resname, "_AddDom_GRM.txt"),row.names=F,col.names=T,quote=F)

plotdata1 = cbind(grmAdd[,4],grmDom[,4])
colors <- densCols(plotdata1)
png(file="Plot1_AddDom_GRM.png", width=6000, height=6000, res=600, type="cairo")
par(mai=c(1.2,1.3,0.5,0.5), mgp=c(3.4,1.2,0))
plot(plotdata1, xlab="Off-diagonal elements of additive GRM", ylab="Off-diagonal elements of dominance GRM", col=colors, cex.lab=2.2, cex.axis=2, cex=1.7, pch=16)
legend("top", legend=c(paste("cor(Add_GRM,Dom_GRM)=",round(cor(plotdata1[,1],plotdata1[,2]),3),sep="")), bty="n", cex=2.2)
dev.off()


######################################################################################################################
### Plot2: Scatter Plot of all traits using phenotypic variation explained by additive effect & dominance effect
######################################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

outdir = "/SAN/mottlab/heterosis/3_HSrats/5VC_gene/HSrats_heart_Trpts"; setwd(outdir)
resname = "HSrats_heart"

raw_dat = read.table(paste0(resname, "_AddDom_Variance.txt"), header=T)
raw_dat[,1] = gsub("vc_", "", raw_dat[,1]); rownames(raw_dat) = raw_dat[,1]

col_dat = cbind("Trait"=raw_dat[,1], "Trait_Class"="Class1", "Color"="black")
write.table(col_dat, file=paste0(resname, "_TrptVC_Colors.txt"), row.names=F, col.names=T, quote=F)
col_dat = read.table(paste0(resname, "_TrptVC_Colors.txt"), header=T); rownames(col_dat) = col_dat[,1]

plotdata2 = raw_dat[,c("Add_Var","Dom_Var")]; rownames(plotdata2) = raw_dat[,1]
plotdata2 = cbind(plotdata2,cols="black"); plotdata2[,3] = col_dat[rownames(plotdata2),3]
plotdata2[,3][is.na(plotdata2[,3])] = "black"

png(file="Plot2_AddDom_Variance.png", width=13000, height=8500, res=600, type="cairo")
par(mai=c(1.3,1.3,0.6,0.6), mgp=c(3.5,1.2,0))
xylim = max(as.numeric(plotdata2[,1]),as.numeric(plotdata2[,2]))*1.1
plot(plotdata2[,1],plotdata2[,2],xlab = "Phenotypic variation explained by additive effect",ylab = "Phenotypic variation explained by dominance effect",col=plotdata2[,3],cex.lab =2.7,cex.axis=2.6,cex=2.2,pch=19,xlim = c(0,xylim),ylim=c(0,xylim))
abline(a=0,b=1)
colors_legend = col_dat[!duplicated(col_dat[,2]),]; colors_legend = colors_legend[order(colors_legend[,2]),]
legend("topleft",inset=.05,title="Traits class",legend=colors_legend[,2],col=colors_legend[,3],pch=rep(19,),cex=1.8)
dev.off()


######################################################################################################################
### Plot3: Boxplot of key traits using phenotypic variation explained by additive effect & dominance effect
######################################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

outdir = "/SAN/mottlab/heterosis/3_HSrats/5VC_gene/HSrats_heart_Trpts"; setwd(outdir)
resname = "HSrats_heart"
thres_vc = 0.8

data = read.table(file=paste0(resname, "_AddDom_Variance.txt"), header=T)
data[,1] = gsub("vc_", "", data[,1]); rownames(data) = data[,1]
data = subset(data, data[,"Add_Var"]>thres_vc | data[,"Dom_Var"]>thres_vc)

traits_plot = data[,1]; plotdata3 = data[traits_plot,]
plotdata3 = cbind(plotdata3[,1:3],Add_col="blue",Add_lty=1,plotdata3[4:5],Dom_col="red",Dom_lty=4)

# Rearrange the individual order #
tmp1 = subset(plotdata3,plotdata3[,"Add_Var"] > plotdata3[,"Dom_Var"])
tmp1 = tmp1[order(-tmp1[,"Add_Var"]),]
tmp2 = subset(plotdata3,plotdata3[,"Add_Var"] < plotdata3[,"Dom_Var"])
tmp2 = tmp2[order(-tmp2[,"Dom_Var"]),]; tmp2 = tmp2[,c(1,6,7,8,9,2,3,4,5)]

# Rearrange the Add_Var+Dom_Var to Big_Var+Small_Var #
colnames(tmp1) = c("Trait","Big_Var","Big_SE","Big_col","Big_lty","Small_Var","Small_SE","Small_col","Small_lty")
colnames(tmp2) = c("Trait","Big_Var","Big_SE","Big_col","Big_lty","Small_Var","Small_SE","Small_col","Small_lty")
plotdata3_order = rbind(tmp1,tmp2)
plotdata3_order[,"Big_Var"] = as.numeric(plotdata3_order[,"Big_Var"]); plotdata3_order[,"Big_SE"] = as.numeric(plotdata3_order[,"Big_SE"])
plotdata3_order[,"Small_Var"] = as.numeric(plotdata3_order[,"Small_Var"]); plotdata3_order[,"Small_SE"] = as.numeric(plotdata3_order[,"Small_SE"])

png(file="Plot3_AddDom_VarSE.png",width=30000,height=8000,res=600,type="cairo")
par(mar = c(8, 6, 2, 0) + 0.1, mgp=c(3.4,1,0))
plotTop <- max(c(plotdata3_order[,"Big_Var"],plotdata3_order[,"Small_Var"])) * 1.3

barCenters <- barplot(height = plotdata3_order[,"Big_Var"],
				  names.arg = rownames(plotdata3_order),
                  beside = true, las = 2,
                  ylim = c(0, plotTop),
                  cex.names = 0.75, xaxt = "n",
                  ylab = "Proportion of variance explained",
                  border = "black", axes = TRUE, col = plotdata3_order[,"Big_col"],
				  cex.lab =2, cex.axis=1.7)
text(x = barCenters, y = -0.02, srt = 45, adj = 1, labels = rownames(plotdata3_order), xpd = TRUE, cex = 1)
segments(barCenters, plotdata3_order[,"Big_Var"], barCenters, plotdata3_order[,"Big_Var"]+plotdata3_order[,"Big_SE"], lwd = 1.5, lty = plotdata3_order[,"Big_lty"])
arrows(barCenters, plotdata3_order[,"Big_Var"], barCenters, plotdata3_order[,"Big_Var"]+plotdata3_order[,"Big_SE"], lwd = 1.5, angle = 90, length = 0.05, lty = plotdata3_order[,"Big_lty"])  	

barCenters <- barplot(height = plotdata3_order[,"Small_Var"],
				  names.arg = rownames(plotdata3_order),
                  beside = true, las = 2,
                  ylim = c(0, plotTop),
                  cex.names = 0.75, xaxt = "n",
                  ylab = "Proportion of variance explained",
                  border = "black", axes = TRUE, col = plotdata3_order[,"Small_col"], add = T,
				  cex.lab =2, cex.axis=1.7)
text(x = barCenters, y = -0.02, srt = 45, adj = 1, labels = rownames(plotdata3_order), xpd = TRUE, cex = 1)
segments(barCenters, plotdata3_order[,"Small_Var"], barCenters, plotdata3_order[,"Small_Var"]+plotdata3_order[,"Small_SE"], lwd = 1.5, lty=plotdata3_order[,"Small_lty"])
arrows(barCenters, plotdata3_order[,"Small_Var"], barCenters, plotdata3_order[,"Small_Var"]+plotdata3_order[,"Small_SE"], lwd = 1.5, angle = 90, length = 0.05, lty=plotdata3_order[,"Small_lty"])

legend("top",bty="n",legend=c(expression(delta^2~"SNP"),expression(h^2~"SNP")),col = c("red","blue"),pch=15,cex=3)  

dev.off()

