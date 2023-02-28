
rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/SAN/mottlab/heterosis/1_F2pigs/3NonAdd_QTL/1_A+D/2_Pvalue"); chr_rm = c(23)
file_all  = dir(); SizeMatrix = length(file_all)
color_sum = read.table("/SAN/mottlab/heterosis/1_F2pigs/1data/F2.cols",header=T,sep="\t")

# png_title = "tODtA_AvsAD_Plot"; pchoose_col = "X.logP.AddDomCode_Add."; model_name = "AvsAD"; stock_name = "F2 Pigs"
# png_title = "tODtA_D_Plot"; pchoose_col = "X.logP.DomCode."; model_name = "D"; stock_name = "F2 Pigs"
# png_title = "tODtA_A_Plot"; pchoose_col = "X.logP.AddCode."; model_name = "A"; stock_name = "F2 Pigs"
png_title = "tODtA_AD_Plot"; pchoose_col = "X.logP.AddDomCode_Null."; model_name = "AD"; stock_name = "F2 Pigs"

##################################################################################################
#### [1] Plot1: Summary and Plotting of the tDom/tAdd distribution of all SNPs on all traits
##################################################################################################

png(file=paste("../",png_title,".png",sep=""),width=16000,height=16000,res=100,type="cairo")
layout(matrix(1:(ceiling(sqrt(SizeMatrix))*ceiling(sqrt(SizeMatrix))),nrow=ceiling(sqrt(SizeMatrix))),widths=rep(1,16),heights=rep(1,16))
par(mai=c(1.4,1.4,1.2,1),mgp=c(3.6,1.3,0))

rownames(color_sum) = color_sum[,1]; color_sum_reorder = color_sum[order(color_sum[,2],color_sum[,1]),]
file_all_reorder = color_sum[gsub(".txt.tar.gz","",gsub("Pvalue_","",file_all)),1] # Reorder the phenotypes by the "Trait_Class" and "Trait" #

tmp_plotdata_all = c()
for(i in 1:length(file_all_reorder)){
try({
	tmp_traitname = file_all_reorder[i]
	tmp_filename = paste("Pvalue_",tmp_traitname,".txt.tar.gz",sep="")
	system(paste("tar zxvf ",tmp_filename,sep=""))
	tmp_data = read.table(file=paste(gsub(".tar.gz","",tmp_filename),sep=""),header=T)
	tmp_data = subset(tmp_data,tmp_data[,1]!=0 | tmp_data[,2]!=0)
	tmp_data = tmp_data[complete.cases(tmp_data),]
	system(paste("rm ",gsub(".tar.gz","",tmp_filename),sep=""))
	
	thrgenome = -log10(0.05/nrow(tmp_data)); thrsuggest = -log10(1/nrow(tmp_data))
	tmp_data_suggest = subset(tmp_data,tmp_data[,pchoose_col] > thrsuggest)
	tmp_data_genome = subset(tmp_data,tmp_data[,pchoose_col] > thrgenome)
	
	colors <- densCols(tmp_data[,c("tDom","tAdd")])
	xylim = max(abs(c(tmp_data[,"tDom"],tmp_data[,"tAdd"])))*1.05
	plot(tmp_data[,c("tDom","tAdd")],col=colors,pch=20,main=tmp_traitname,xlab="t(OverDom)",ylab="t(Add)",xlim=c(-xylim,xylim),ylim=c(-xylim,xylim),cex=2.5,cex.axis=2.3,cex.lab=2.5,cex.main=4.5)
	points(tmp_data_suggest[,c("tDom","tAdd")],col="red",pch=20,cex=2.7)
	points(tmp_data_genome[,c("tDom","tAdd")],col="red4",pch=20,cex=2.7)
	legend("topleft",legend=c(paste("cor(X,Y) = ",round(cor(tmp_data[,"tDom"],tmp_data[,"tAdd"]),3),sep="")),bty="n",cex=3.5)
	
	try({
		length2 = xylim*1.5
		segments(-length2,-length2/1.2,length2,length2/1.2,lty=2,lwd=1.3,col="gray17"); segments(-length2*0.8,-length2,length2*0.8,length2,lty=2,lwd=1.3,col="gray17"); segments(-length2*0.2,-length2,length2*0.2,length2,lty=2,lwd=1.3,col="gray17")
		segments(-length2,length2/1.2,length2,-length2/1.2,lty=2,lwd=1.3,col="gray17"); segments(-length2*0.8,length2,length2*0.8,-length2,lty=2,lwd=1.3,col="gray17"); segments(-length2*0.2,length2,length2*0.2,-length2,lty=2,lwd=1.3,col="gray17")
		arrows(-length2*1.1,0,length2*1.2,0,lwd=2.2,code=2,angle=17); arrows(0,-length2*1.1,0,length2*1.2,lwd=2.2,code=2,angle=17)
	})
	
	if(nrow(tmp_data_suggest) > 0){
		tmp_plotdata = cbind(trait=tmp_traitname,tmp_data_suggest)
		tmp_plotdata_all = rbind(tmp_plotdata_all,tmp_plotdata) 
	}

	print(i)
})
}
dev.off()


##################################################################################################
#### [2] Summary the peak QTLs of each trait and resort by the trait_class
##################################################################################################

### Part1: Get the significant peak on different chrs of each trait ###
SigTraits = unique(tmp_plotdata_all[,"trait"])
SigTraits_Peak = c()
for(i in 1:length(SigTraits)){
	TmpData = subset(tmp_plotdata_all,tmp_plotdata_all[,"trait"]==SigTraits[i])
	TmpPeak_all = c()
	for(j in 1:length(unique(TmpData[,"chr"]))){
		TmpData_Tmpchr = subset(TmpData,TmpData[,"chr"]==unique(TmpData[,"chr"])[j])
		TmpPeak = subset(TmpData_Tmpchr,TmpData_Tmpchr[,pchoose_col]==max(TmpData_Tmpchr[,pchoose_col]))[1,]
		TmpPeak_all = rbind(TmpPeak_all,TmpPeak)
	}
	SigTraits_Peak = rbind(SigTraits_Peak,TmpPeak_all)
}
SigTraits_Peak = subset(SigTraits_Peak,SigTraits_Peak[,"chr"]!=chr_rm)

### Part2: Resort the trait by the trait_class ###
for(j in 1:nrow(SigTraits_Peak)){ SigTraits_Peak[j,"Trait_Class"] = color_sum[SigTraits_Peak[j,"trait"],"Trait_Class"] }
SigTraits_Peak = SigTraits_Peak[order(SigTraits_Peak[,"Trait_Class"]),]
save(SigTraits_Peak,color_sum,file=paste("../",png_title,"_Sig.RData",sep=""))


##################################################################################################
#### [3] Plot2: Significant QTLs of each traits classified by tDom/tAdd (Circle distribution)
##################################################################################################

try({

	png(file=paste("../",png_title,"_Sig1.png",sep=""),width=7000*1.5,height=7500,res=600,type="cairo")
	layout(matrix(1:2,nrow=1),widths=c(1,0.5),heights=1)
	par(mai=c(1.4,1.4,1.2,0.1),mgp=c(3.6,1.3,0))

	colors <- color_sum[SigTraits_Peak[,"trait"],"Color"]; xylim = max(abs(c(SigTraits_Peak[,"tDom"],SigTraits_Peak[,"tAdd"])))*1.05
	plot(0,0,main=paste0(stock_name,": Significant QTLs detected by ",model_name," model"),xlab="t(Dom)",ylab="t(Add)",xlim=c(-xylim*1.15,xylim*1.15),ylim=c(-xylim*1.15,xylim*1.15),cex=0.01,cex.axis=2.3,cex.lab=2.7,cex.main=2.3,bty='n')

	try({
		length2 = xylim; col_add = "grey23"; col_partial = "grey43"; col_dom = "grey63"; col_over = "grey83"

		polygon(c(-length2*0.2,length2*0.2,0),c(length2,length2,0), xpd=T, col=col_add, border=F)
		polygon(c(-length2*0.2,length2*0.2,0),c(-length2,-length2,0), xpd=T, col=col_add, border=F)

		polygon(c(-length2*0.8,-length2*0.2,0),c(length2,length2,0), xpd=T, col=col_partial, border=F)
		polygon(c(length2*0.8,length2*0.2,0),c(length2,length2,0), xpd=T, col=col_partial, border=F)
		polygon(c(length2*0.8,length2*0.2,0),c(-length2,-length2,0), xpd=T, col=col_partial, border=F)
		polygon(c(-length2*0.8,-length2*0.2,0),c(-length2,-length2,0), xpd=T, col=col_partial, border=F)

		polygon(c(-length2,-length2,0,-length2*0.8),c(length2,length2*0.83,0,length2), xpd=T, col=col_dom, border=F)
		polygon(c(length2,length2,0,length2*0.8),c(length2,length2*0.83,0,length2), xpd=T, col=col_dom, border=F)
		polygon(c(length2,length2,0,length2*0.8),c(-length2,-length2*0.83,0,-length2), xpd=T, col=col_dom, border=F)
		polygon(c(-length2,-length2,0,-length2*0.8),c(-length2,-length2*0.83,0,-length2), xpd=T, col=col_dom, border=F)

		polygon(c(-length2,-length2,0),c(-length2*0.83,length2*0.83,0), xpd=T, col=col_over, border=F)
		polygon(c(length2,length2,0),c(-length2*0.83,length2*0.83,0), xpd=T, col=col_over, border=F)

		segments(-length2,-length2/1.2,length2,length2/1.2,lty=2,lwd=1.3,col="gray17"); segments(-length2*0.8,-length2,length2*0.8,length2,lty=2,lwd=1.3,col="gray17"); segments(-length2*0.2,-length2,length2*0.2,length2,lty=2,lwd=1.3,col="gray17")
		segments(-length2,length2/1.2,length2,-length2/1.2,lty=2,lwd=1.3,col="gray17"); segments(-length2*0.8,length2,length2*0.8,-length2,lty=2,lwd=1.3,col="gray17"); segments(-length2*0.2,length2,length2*0.2,-length2,lty=2,lwd=1.3,col="gray17")
		#segments(-length2,-length2,length2,length2,lty=1,lwd=1.5); segments(-length2,length2,length2,-length2,lty=1,lwd=1.5)
		
		arrows(-length2*1.1,0,length2*1.2,0,lwd=2.2,code=2,angle=17)
		arrows(0,-length2*1.1,0,length2*1.2,lwd=2.2,code=2,angle=17)
		
		legend("topright",c("Additive (0<|tDom/tAdd|<0.2)","Partial-Dominant (0.2<|tDom/tAdd|<0.8)","Dominant (0.8<|tDom/tAdd|<1.2)","Over-Dominant (1.2<|tDom/tAdd|)"),pch=c(15,15,15,15), col=c(col_add,col_partial,col_dom,col_over),cex=1.25)
		
	})
	 
	points(SigTraits_Peak[,c("tDom","tAdd")],col=colors,pch=20,cex=2.5)

	par(mai=c(1,0.01,0.5,0.1),mgp=c(3,2,0),tck=0)
	plot(x=0,y=0,pch=c("."),cex.axis=0.001,cex.lab=0.001,axes=F)
	colors_legend = color_sum[SigTraits_Peak[,"trait"],][!duplicated(color_sum[SigTraits_Peak[,"trait"],2]),]
	colors_legend = colors_legend[order(colors_legend[,2]),]

	legend("center",legend=colors_legend[,2],col=colors_legend[,3],pch=16,bty="n",cex=1)

	dev.off()

})


#####################################################################################################################
#### [4] Plot3: Significant QTLs of each traits classified by tDom/tAdd (Transverse distribution ordered by tDtA)
#####################################################################################################################

try({

	png(file=paste("../",png_title,"_Sig2.png",sep=""),width=18000,height=9000,res=600,type="cairo")
	layout(matrix(1:2,nrow=1),widths=c(1,0.2),heights=1)
	par(mar=c(14,7,7,0),mgp=c(2.7,0.8,0),cex.main=3)

	SigTraits_Peak_order = cbind("order"=log2(abs(SigTraits_Peak[,"tDom"]/SigTraits_Peak[,"tAdd"])), SigTraits_Peak)
	SigTraits_Peak_order = SigTraits_Peak_order[order(SigTraits_Peak_order[,"Trait_Class"],-SigTraits_Peak_order[,"order"]),-1]; SigTraits_Peak = SigTraits_Peak_order
	# SigTraits_Peak_order = SigTraits_Peak_order[order(SigTraits_Peak_order[,"Trait_Class"],-SigTraits_Peak_order[,pchoose_col]),-1]; SigTraits_Peak = SigTraits_Peak_order

	### Part0: Calculate the num and ratio of D-QTLs & OD-QTLs ###
	ratio_value = cbind(SigTraits_Peak[,c("Trait_Class","trait")], "tDtA"=abs(SigTraits_Peak[,"tDom"]/SigTraits_Peak[,"tAdd"]), "qtl_type"="A", "plot_order"=1:nrow(SigTraits_Peak), "Color"=color_sum[SigTraits_Peak[,"trait"],"Color"])
	ratio_value[(ratio_value[,"tDtA"]>0.2 & ratio_value[,"tDtA"]<=0.8), "qtl_type"] = "PD"
	ratio_value[(ratio_value[,"tDtA"]>0.8 & ratio_value[,"tDtA"]<=1.2), "qtl_type"] = "D"
	ratio_value[(ratio_value[,"tDtA"]>1.2), "qtl_type"] = "OD"
	traitclass = unique(ratio_value[,1]); ratio_plot = c()
	for(i in 1:length(traitclass)){ 
		tmp1 = ratio_value[ratio_value[,1]==traitclass[i],]
		tmp2 = tmp1[(tmp1[,"qtl_type"] %in% c("D", "OD")),]
		ratio_plot_tmp = c(traitclass[i], mean(tmp1[,"plot_order"]), nrow(tmp2), nrow(tmp2)/nrow(tmp1), unique(tmp1[,"Color"]))
		ratio_plot = rbind(ratio_plot, ratio_plot_tmp)
	}

	### Part1: The main plot of tDom/tAdd ###
	col_add = "grey23"; col_partial = "grey43"; col_dom = "grey63"; col_over = "grey83"
	plot_x = 1:nrow(SigTraits_Peak)
	plot_y = abs(SigTraits_Peak[,"tDom"]/SigTraits_Peak[,"tAdd"])
	plot_cex = SigTraits_Peak[,"X.logP.AddDomCode_Null."]; zoomout = 0.6
	plot_pch = rep(20, nrow(SigTraits_Peak))
	plot_pch[plot_cex>10] = 18 # 将-logP值大于10的点，用菱形绘制 #
	plot_cex[plot_cex>10] = 8 # 将-logP值大于10的点，大小全部归成10 #
	plot_col = color_sum[SigTraits_Peak[,"trait"],"Color"] # 依据color_sum来给每个性状赋色 #
	plot_ylim = 100; axis_pos = -7.2

	plot(plot_x, log2(plot_y), cex=plot_cex*zoomout, col=plot_col, pch=plot_pch, xaxt="n", xlab="", ylab="|tDom/tAdd| (log2 transformed)", main=paste0(stock_name,": Significant QTLs detected by ",model_name," model"), cex.lab=2.6, cex.axis=2, cex.main=3.7)
	polygon(c(-nrow(SigTraits_Peak),nrow(SigTraits_Peak)*2,nrow(SigTraits_Peak)*2,-nrow(SigTraits_Peak)),log2(c(min(plot_y)/2,min(plot_y)/2,0.2,0.2)),col=col_add,xpd=F,border=F)
	polygon(c(-nrow(SigTraits_Peak),nrow(SigTraits_Peak)*2,nrow(SigTraits_Peak)*2,-nrow(SigTraits_Peak)),log2(c(0.2,0.2,0.8,0.8)),col=col_partial,xpd=F,border=F)
	polygon(c(-nrow(SigTraits_Peak),nrow(SigTraits_Peak)*2,nrow(SigTraits_Peak)*2,-nrow(SigTraits_Peak)),log2(c(0.8,0.8,1.2,1.2)),col=col_dom,xpd=F,border=F)
	polygon(c(-nrow(SigTraits_Peak),nrow(SigTraits_Peak)*2,nrow(SigTraits_Peak)*2,-nrow(SigTraits_Peak)),log2(c(1.2,1.2,plot_ylim*2,plot_ylim*2)),col=col_over,xpd=F,border=F)
	points(plot_x, log2(plot_y), cex=plot_cex*zoomout, col=plot_col, pch=plot_pch)
	# text(x = plot_x, y = axis_pos, srt = 45, adj = 1, labels = gsub(".resid","",SigTraits_Peak[,"trait"]), xpd = TRUE, cex = 0.5)
	text(x = plot_x[seq(1,length(plot_x),2)], y = axis_pos, srt = 90, adj = 1, labels = gsub(".resid","",SigTraits_Peak[,"trait"])[seq(1,length(plot_x),2)], xpd = TRUE, cex = 0.8)
	text(x = plot_x[seq(2,length(plot_x),2)], y = axis_pos-1.3, srt = 90, adj = 1, labels = gsub(".resid","",SigTraits_Peak[,"trait"])[seq(2,length(plot_x),2)], xpd = TRUE, cex = 0.8)
	text(x = as.numeric(ratio_plot[,2]), y = axis_pos+0.5, adj = 0.5, labels = paste0(ratio_plot[,3], " (", signif(as.numeric(ratio_plot[,4])*100,3), "%)"), col = ratio_plot[,5], xpd = TRUE, cex = 1.5)
	
	### Part2: The legend plot of trait colors, p-value size, tDom/tAdd region size ###
	par(mai=c(0.1,0.01,0.5,0.1),mgp=c(3,2,0),tck=0)
	plot(x=0,y=0,pch=c("."),cex.axis=0.001,cex.lab=0.001,axes=F)
	colors_legend = color_sum[SigTraits_Peak[,"trait"],][!duplicated(color_sum[SigTraits_Peak[,"trait"],2]),]
	colors_legend = colors_legend[order(colors_legend[,2]),]
	legend("center",legend=colors_legend[,2],col=colors_legend[,3],pch=20,cex=1,bty="o")

	legend("top",c("Over-Dominant QTL","(1.2<|tDom/tAdd|)","Dominant QTL","(0.8<|tDom/tAdd|<1.2)","Partial-Dominant QTL","(0.2<|tDom/tAdd|<0.8)","Additive QTL","(0<|tDom/tAdd|<0.2)"),pch=c(15,15,15,15,15,15,15,15), col=c(col_over,"white",col_dom,"white",col_partial,"white",col_add,"white"),cex=1,bty="o")

	legend("bottom",legend=c("  -logP = 4","","  -logP = 6","","  -logP = 8","","  -logP > 10",""),col=c("black","white","black","white","black","white","black","white"),pch=c(20,20,20,20,20,20,18,20),pt.cex=c(4*zoomout,0.01,6*zoomout,0.01,8*zoomout,0.01,8*zoomout,0.01),cex=1.7,bty="o")

	dev.off()

})


############################################################################################################################
#### [5] Plot4: Significant QTLs of each traits classified by tDom/tAdd (Transverse distribution ordered by pchoose_col)
############################################################################################################################

try({

	png(file=paste("../",png_title,"_Sig3.png",sep=""),width=18000,height=9000,res=600,type="cairo")
	layout(matrix(1:2,nrow=1),widths=c(1,0.2),heights=1)
	par(mar=c(14,7,7,0),mgp=c(2.7,0.8,0),cex.main=3)

	SigTraits_Peak_order = cbind("order"=log2(abs(SigTraits_Peak[,"tDom"]/SigTraits_Peak[,"tAdd"])), SigTraits_Peak)
	# SigTraits_Peak_order = SigTraits_Peak_order[order(SigTraits_Peak_order[,"Trait_Class"],-SigTraits_Peak_order[,"order"]),-1]; SigTraits_Peak = SigTraits_Peak_order
	SigTraits_Peak_order = SigTraits_Peak_order[order(SigTraits_Peak_order[,"Trait_Class"],-SigTraits_Peak_order[,pchoose_col]),-1]; SigTraits_Peak = SigTraits_Peak_order

	### Part0: Calculate the num and ratio of D-QTLs & OD-QTLs ###
	ratio_value = cbind(SigTraits_Peak[,c("Trait_Class","trait")], "tDtA"=abs(SigTraits_Peak[,"tDom"]/SigTraits_Peak[,"tAdd"]), "qtl_type"="A", "plot_order"=1:nrow(SigTraits_Peak), "Color"=color_sum[SigTraits_Peak[,"trait"],"Color"])
	ratio_value[(ratio_value[,"tDtA"]>0.2 & ratio_value[,"tDtA"]<=0.8), "qtl_type"] = "PD"
	ratio_value[(ratio_value[,"tDtA"]>0.8 & ratio_value[,"tDtA"]<=1.2), "qtl_type"] = "D"
	ratio_value[(ratio_value[,"tDtA"]>1.2), "qtl_type"] = "OD"
	traitclass = unique(ratio_value[,1]); ratio_plot = c()
	for(i in 1:length(traitclass)){ 
		tmp1 = ratio_value[ratio_value[,1]==traitclass[i],]
		tmp2 = tmp1[(tmp1[,"qtl_type"] %in% c("D", "OD")),]
		ratio_plot_tmp = c(traitclass[i], mean(tmp1[,"plot_order"]), nrow(tmp2), nrow(tmp2)/nrow(tmp1), unique(tmp1[,"Color"]))
		ratio_plot = rbind(ratio_plot, ratio_plot_tmp)
	}

	### Part1: The main plot of tDom/tAdd ###
	col_add = "grey23"; col_partial = "grey43"; col_dom = "grey63"; col_over = "grey83"
	plot_x = 1:nrow(SigTraits_Peak)
	plot_y = abs(SigTraits_Peak[,"tDom"]/SigTraits_Peak[,"tAdd"])
	plot_cex = SigTraits_Peak[,"X.logP.AddDomCode_Null."]; zoomout = 0.6
	plot_pch = rep(20, nrow(SigTraits_Peak))
	plot_pch[plot_cex>10] = 18 # 将-logP值大于10的点，用菱形绘制 #
	plot_cex[plot_cex>10] = 8 # 将-logP值大于10的点，大小全部归成10 #
	plot_col = color_sum[SigTraits_Peak[,"trait"],"Color"] # 依据color_sum来给每个性状赋色 #
	plot_ylim = 100; axis_pos = -7.2

	plot(plot_x, log2(plot_y), cex=plot_cex*zoomout, col=plot_col, pch=plot_pch, xaxt="n", xlab="", ylab="|tDom/tAdd| (log2 transformed)", main=paste0(stock_name,": Significant QTLs detected by ",model_name," model"), cex.lab=2.6, cex.axis=2, cex.main=3.7)
	polygon(c(-nrow(SigTraits_Peak),nrow(SigTraits_Peak)*2,nrow(SigTraits_Peak)*2,-nrow(SigTraits_Peak)),log2(c(min(plot_y)/2,min(plot_y)/2,0.2,0.2)),col=col_add,xpd=F,border=F)
	polygon(c(-nrow(SigTraits_Peak),nrow(SigTraits_Peak)*2,nrow(SigTraits_Peak)*2,-nrow(SigTraits_Peak)),log2(c(0.2,0.2,0.8,0.8)),col=col_partial,xpd=F,border=F)
	polygon(c(-nrow(SigTraits_Peak),nrow(SigTraits_Peak)*2,nrow(SigTraits_Peak)*2,-nrow(SigTraits_Peak)),log2(c(0.8,0.8,1.2,1.2)),col=col_dom,xpd=F,border=F)
	polygon(c(-nrow(SigTraits_Peak),nrow(SigTraits_Peak)*2,nrow(SigTraits_Peak)*2,-nrow(SigTraits_Peak)),log2(c(1.2,1.2,plot_ylim*2,plot_ylim*2)),col=col_over,xpd=F,border=F)
	points(plot_x, log2(plot_y), cex=plot_cex*zoomout, col=plot_col, pch=plot_pch)
	# text(x = plot_x, y = axis_pos, srt = 45, adj = 1, labels = gsub(".resid","",SigTraits_Peak[,"trait"]), xpd = TRUE, cex = 0.5)
	text(x = plot_x[seq(1,length(plot_x),2)], y = axis_pos, srt = 90, adj = 1, labels = gsub(".resid","",SigTraits_Peak[,"trait"])[seq(1,length(plot_x),2)], xpd = TRUE, cex = 0.8)
	text(x = plot_x[seq(2,length(plot_x),2)], y = axis_pos-1.3, srt = 90, adj = 1, labels = gsub(".resid","",SigTraits_Peak[,"trait"])[seq(2,length(plot_x),2)], xpd = TRUE, cex = 0.8)
	text(x = as.numeric(ratio_plot[,2]), y = axis_pos+0.5, adj = 0.5, labels = paste0(ratio_plot[,3], " (", signif(as.numeric(ratio_plot[,4])*100,3), "%)"), col = ratio_plot[,5], xpd = TRUE, cex = 1.5)
	
	### Part2: The legend plot of trait colors, p-value size, tDom/tAdd region size ###
	par(mai=c(0.1,0.01,0.5,0.1),mgp=c(3,2,0),tck=0)
	plot(x=0,y=0,pch=c("."),cex.axis=0.001,cex.lab=0.001,axes=F)
	colors_legend = color_sum[SigTraits_Peak[,"trait"],][!duplicated(color_sum[SigTraits_Peak[,"trait"],2]),]
	colors_legend = colors_legend[order(colors_legend[,2]),]
	legend("center",legend=colors_legend[,2],col=colors_legend[,3],pch=20,cex=1,bty="o")

	legend("top",c("Over-Dominant QTL","(1.2<|tDom/tAdd|)","Dominant QTL","(0.8<|tDom/tAdd|<1.2)","Partial-Dominant QTL","(0.2<|tDom/tAdd|<0.8)","Additive QTL","(0<|tDom/tAdd|<0.2)"),pch=c(15,15,15,15,15,15,15,15), col=c(col_over,"white",col_dom,"white",col_partial,"white",col_add,"white"),cex=1,bty="o")

	legend("bottom",legend=c("  -logP = 4","","  -logP = 6","","  -logP = 8","","  -logP > 10",""),col=c("black","white","black","white","black","white","black","white"),pch=c(20,20,20,20,20,20,18,20),pt.cex=c(4*zoomout,0.01,6*zoomout,0.01,8*zoomout,0.01,8*zoomout,0.01),cex=1.7,bty="o")

	dev.off()

})
