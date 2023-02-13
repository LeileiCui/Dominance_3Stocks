
### STEP3. Plot ###

######################################################################################################################
### Plot1: Scatter plot of all individuals-pairs using off-diagonal elements of additive GRM & dominance GRM
######################################################################################################################
rm(list=ls());options(stringsAsFactors=FALSE)
outdir = "/SAN/mottlab/heterosis/4_HSmice/5VC/AllInds/"; setwd(outdir)

grmAdd = read.table(file="HSmice_Final_add.grm.gz"); grmDom = read.table(file="HSmice_Final_dom.d.grm.gz")
grmAdd = grmAdd[grmAdd[,1]!=grmAdd[,2],]; grmDom = grmDom[grmDom[,1]!=grmDom[,2],]
write.table(cbind(indpair=paste(grmAdd[,1],"_",grmAdd[,2],sep=""),grmAdd=grmAdd[,4],grmDom=grmDom[,4]),file=paste(outdir,"Plot1_AddDom_GRM.txt",sep=""),row.names=F,col.names=T,quote=F)

plotdata1 = cbind(grmAdd[,4],grmDom[,4])
colors <- densCols(plotdata1)
png(file=paste(outdir,"Plot1_AddDom_GRM.png",sep=""),width=6000,height=6000,res=600,type="cairo")
par(mai=c(1.2,1.3,0.5,0.5),mgp=c(3.4,1.2,0))
plot(plotdata1,xlab="Off-diagonal elements of additive GRM",ylab="Off-diagonal elements of dominance GRM",col=colors,cex.lab=2.2,cex.axis=2,cex=1.7,pch=16)
legend("top",legend=c(paste("cor(Add_GRM,Dom_GRM)=",round(cor(plotdata1[,1],plotdata1[,2]),3),sep="")),bty="n",cex=2.2)
dev.off()

######################################################################################################################
### Plot2: Scatter Plot of all traits using phenotypic variation explained by additive effect & dominance effect
######################################################################################################################
rm(list=ls());options(stringsAsFactors=FALSE)
outdir = "/SAN/mottlab/heterosis/4_HSmice/5VC/AllInds/"; setwd(outdir)
outdir_VC = "/SAN/mottlab/heterosis/4_HSmice/5VC/AllInds/add_dom_VC"
hsq_list = list.files(path="/SAN/mottlab/heterosis/4_HSmice/5VC/AllInds/add_dom_VC", pattern=".hsq")

tmp_result_all = c()
for(i in 1:length(hsq_list)){
	tmp_hsq = read.table(paste(outdir_VC,"/",hsq_list[i],sep=""),header=T,fill=T)
	rownames(tmp_hsq) = tmp_hsq[,1]
	tmp_result = c(hsq_list[i],as.numeric(as.character(tmp_hsq["V(G1)/Vp","Variance"])),as.numeric(as.character(tmp_hsq["V(G1)/Vp","SE"])),as.numeric(as.character(tmp_hsq["V(G2)/Vp","Variance"])),as.numeric(as.character(tmp_hsq["V(G2)/Vp","SE"])))
	tmp_result_all = rbind(tmp_result_all,tmp_result)
	print(paste(i," ",hsq_list[i],sep=""))
}
tmp_result_all[,1] = gsub(".hsq","",tmp_result_all[,1])
colnames(tmp_result_all) = c("Trait","Add_Var","Add_SE","Dom_Var","Dom_SE")
#write.table(tmp_result_all,file=paste(outdir,"Plot2_AddDom_Variance.txt",sep=""),row.names=F,col.names=T,quote=F)

color = read.table("/SAN/mottlab/heterosis/4_HSmice/5VC/data/Data_clean_colors.txt",header=T,stringsAsFactors=F,sep="\t"); rownames(color) = color[,1]
plotdata2 = tmp_result_all[,c("Add_Var","Dom_Var")]; rownames(plotdata2) = tmp_result_all[,1]
plotdata2 = cbind(plotdata2,cols="black"); plotdata2[,3] = color[rownames(plotdata2),3]
plotdata2[,3][is.na(plotdata2[,3])] = "black"

png(file=paste(outdir,"Plot2_AddDom_Variance.png",sep=""),width=13000,height=8500,res=600,type="cairo")
par(mai=c(1.3,1.3,0.6,0.6),mgp=c(3.5,1.2,0))
xylim = max(as.numeric(plotdata2[,1]),as.numeric(plotdata2[,2]))*1.1
plot(plotdata2[,1],plotdata2[,2],xlab = "Phenotypic variation explained by additive effect",ylab = "Phenotypic variation explained by dominance effect",col=plotdata2[,3],cex.lab =2.7,cex.axis=2.6,cex=2.2,pch=19,xlim = c(0,xylim),ylim=c(0,xylim))
abline(a=0,b=1)
colors_legend = color[!duplicated(color[,2]),]; colors_legend = colors_legend[order(colors_legend[,2]),]
legend("topright",inset=.05,title="Traits class",legend=colors_legend[,2],col=colors_legend[,3],pch=rep(19,),cex=2)
dev.off()


	rm(list=ls());options(stringsAsFactors=FALSE)
	outdir = "/SAN/mottlab/heterosis/4_HSmice/5VC/AllInds/"; setwd(outdir)

	tmp_result_all = read.table("Plot2_AddDom_Variance.txt", header=T)
	color = read.table("/SAN/mottlab/heterosis/4_HSmice/5VC/data/Data_clean_colors.txt",header=T,stringsAsFactors=F,sep="\t"); rownames(color) = color[,1]
	plotdata2 = tmp_result_all[,c("Add_Var","Dom_Var")]; rownames(plotdata2) = tmp_result_all[,1]
	plotdata2 = cbind(plotdata2,cols="black"); plotdata2[,3] = color[rownames(plotdata2),3]
	plotdata2[,3][is.na(plotdata2[,3])] = "black"

	png(file=paste(outdir,"Plot2_AddDom_Variance_new1.png",sep=""),width=9000,height=9000,res=600,type="cairo")
	nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(5,1), c(1,5), TRUE)

	par(mai=c(1.3,1.3,0.6,0.6),mgp=c(3.5,1.2,0))
	xylim = max(as.numeric(plotdata2[,1]),as.numeric(plotdata2[,2]))*1.1
	plot(plotdata2[,1],plotdata2[,2],xlab = "Phenotypic variation explained by additive effect",ylab = "Phenotypic variation explained by dominance effect",col=plotdata2[,3],cex.lab =2.7,cex.axis=2.6,cex=2.2,pch=19,xlim = c(0,xylim),ylim=c(0,xylim))
	abline(a=0,b=1)

	### Prepare the ratio to zoom out the yhist###
	yhist_ratio = max(as.numeric(plotdata2[,1]))/max(as.numeric(plotdata2[,2]))

	### Prepare dataset for the histogram plot ###
	x <- as.numeric(plotdata2[,1]); y <- as.numeric(plotdata2[,2])
	x_min = min(x)*0.9; x_max = max(x)*1.1 # As x is the distance, so all the value will be positive #
	y_min = min(y)*0.9; y_max = max(y)*1.1 # As y is the logP, so all the value will be positive #
	xrange <- c(x_min,x_max); yrange <- c(y_min,y_max)
	xhist <- hist(x, breaks=seq(x_min,x_max,(abs(x_max)+abs(x_min))/40), plot=FALSE)
	yhist <- hist(y, breaks=seq(y_min,y_max,(abs(y_max)+abs(y_min))/40), plot=FALSE)
	top <- max(c(xhist$counts, yhist$counts))

	par(mar=c(0,7,1,1), mgp=c(4,1.3,0), lwd=1)
	barplot(xhist$counts, axes=T, ylim=c(0, top), space=0, cex.axis=2.3, col="#0000FFB3")

	par(mar=c(0,7,1,1), mgp=c(4,1.3,0), lwd=3)
	barplot(xhist$counts, axes=T, ylim=c(0, top), space=0, cex.axis=2.3, col="white", angle=45, density=7, add=T)

	par(mar=c(7,0,1,1),mgp=c(4,1.3,0), lwd=1)
	yhist_ylim = yhist_ratio*length(yhist$counts)
	barplot(yhist$counts, axes=T, xlim=c(0, top), ylim=c(0, yhist_ylim), space=0, horiz=TRUE, cex.axis=2.3, col="#FF0000B3")

	par(mar=c(7,0,1,1),mgp=c(4,1.3,0), lwd=3)
	yhist_ylim = yhist_ratio*length(yhist$counts)
	barplot(yhist$counts, axes=T, xlim=c(0, top), ylim=c(0, yhist_ylim), space=0, horiz=TRUE, cex.axis=2.3, col="white", angle=45, density=7, add=T)

	dev.off()


######################################################################################################################
### Plot3: Boxplot of key traits using phenotypic variation explained by additive effect & dominance effect
######################################################################################################################
rm(list=ls());options(stringsAsFactors=FALSE)
outdir = "/SAN/mottlab/heterosis/4_HSmice/5VC/AllInds/"
data = read.table(file="Plot2_AddDom_Variance.txt",header=T); rownames(data) = data[,1]
traits_plot = data[,1]
plotdata3 = data[traits_plot,]; plotdata3 = cbind(plotdata3[,1:3],Add_col="blue",Add_lty=1,plotdata3[4:5],Dom_col="red",Dom_lty=4)

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

png(file=paste(outdir,"Plot3_AddDom_VarSE.png",sep=""),width=30000,height=8000,res=600,type="cairo")
par(mar = c(12, 6, 2, 0) + 0.1, mgp=c(3.4,1,0))
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



