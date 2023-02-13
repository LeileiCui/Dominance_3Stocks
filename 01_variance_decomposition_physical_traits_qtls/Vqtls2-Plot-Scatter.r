
rm(list=ls());options(stringsAsFactors=FALSE)

dir=c("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/2-全基因组方差剖分-2散点图QTL")

for(i in 1:3){

	setwd(paste0(dir, c("/1-F2pigs", "/2-HSrats", "/3-HSmice")[i]))
	files=c("1-Pigs-Vqtl.txt","2-Rats-Vqtl.txt","3-Mice-Vqtl.txt")
	filesCols=c("1-Pigs-Cols.txt","2-Rats-Cols.txt","3-Mice-Cols.txt")
	filenames=c("Pigs","Rats","Mice")
	
	plotdata = read.table(files[i], header=T); plotdata = plotdata[,-c(4,5)]
	plotdata = subset(plotdata, plotdata[,2]<1 & plotdata[,3]<1)
	color = read.table(filesCols[i], header=T, sep="\t"); rownames(color) = color[,1]
	plotdata = cbind(plotdata,cols="black"); plotdata[,8] = color[plotdata[,1],3]; plotdata[,8][is.na(plotdata[,8])] = "black"

	# png(file=paste0(filenames[i],".png"),width=13000,height=8500,res=600,type="cairo")
	# par(mai=c(1.3,1.3,0.6,0.6),mgp=c(3.5,1.2,0))
	# xylim = max(as.numeric(plotdata[,2]),as.numeric(plotdata[,3]))*1.1
	# plot(plotdata[,2],plotdata[,3],xlab = "Phenotypic variation explained by additive effect",ylab = "Phenotypic variation explained by dominance effect",
	# 	col=plotdata[,8],cex.lab =2.7,cex.axis=2.6,cex=2.2,pch=19,xlim = c(0,xylim),ylim=c(0,xylim))
	# abline(a=0,b=1)
	# color = color[(color[,1] %in% unique(plotdata[,1])),]
	# colors_legend = color[!duplicated(color[,2]),]; colors_legend = colors_legend[order(colors_legend[,2]),]
	# legend("topright",inset=.05,title="Traits class",legend=colors_legend[,2],col=colors_legend[,3],pch=rep(19,),cex=1.7)
	# dev.off()

		plotdata2 = plotdata[,c("Vqtla_Vp", "Vqtld_Vp", "cols")]
		png(file=paste0(filenames[i],"_new1.png"),width=9000,height=9000,res=600,type="cairo")
		nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(5,1), c(1,5), TRUE)

		par(mai=c(1.3,1.3,0.6,0.6),mgp=c(3.5,1.2,0))
		xylim = max(as.numeric(plotdata2[,1]),as.numeric(plotdata2[,2]))*1.1
		plot(plotdata2[,1],plotdata2[,2],xlab = "Phenotypic variation explained by additive effect of QTLs",ylab = "Phenotypic variation explained by dominance effect of QTLs",col=plotdata2[,3],cex.lab =2.7,cex.axis=2.6,cex=2.2,pch=19,xlim = c(0,xylim),ylim=c(0,xylim))
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


	print(i)

}
