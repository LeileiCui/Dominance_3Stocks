
rm(list=ls());options(stringsAsFactors=FALSE)

dir=c("/Users/leileicui/Desktop/1_项目资料/2_NonAdditive/结果汇总/Heterosis_Results_2【初步整理结果】/2-全基因组方差剖分-2")

for(i in 1:3){

	setwd(paste0(dir, c("/1-F2pigs", "/2-HSrats", "/3-HSmice")[i]))
	files=c("1-Pigs-Vqtl.txt","2-Rats-Vqtl.txt","3-Mice-Vqtl.txt"); filenames=c("Pigs","Rats","Mice")
	plotdata = read.table(files[i], header=T)[,2:3]
	plotdata = subset(plotdata, plotdata[,1]<1 & plotdata[,2]<1)
	
	png(paste0(filenames[i],"_hist.png"),width=1500,height=1500,res=200,type="cairo")
	par(mgp=c(2.4,0.6,0))

	h1 = hist(plotdata[,1], breaks=seq(0,max(plotdata)+0.05,0.01), plot=FALSE)
	h2 = hist(plotdata[,2], breaks=seq(0,max(plotdata)+0.05,0.01), plot=FALSE)
	h2$counts = - h2$counts
	hmax = max(h1$counts)
	hmin = min(h2$counts)
	X = c(h1$breaks, h2$breaks)
	xmax = max(X)
	xmin = min(X)

	plot(h1, ylim=c(hmin, hmax), col="blue", xlim=c(xmin, xmax), cex.main=1.7, cex.lab=1.7,cex.axis=1.5, 
		main=paste0("Histogram of Va and Vd of all QTLs (",filenames[i],")"), xlab="Phenotypic variance explained by one QTL", ylab="Count")
	lines(h2, col="red")
	legend("topright", legend=c("Additive", "Dominance"), col=c("blue", "red"), border=F, bty="n", lty=1, lwd=7, cex=1.3)
	box()
	dev.off()
	
	print(i)

}
