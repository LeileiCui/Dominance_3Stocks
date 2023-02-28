
rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/home/ucbtlcu/heterosis_mnt/1_F2pigs/3NonAdd_QTL/1_A+D/2_Pvalue"); chr_rm = c(23)
color_sum = read.table("/home/ucbtlcu/heterosis_mnt/1_F2pigs/1data/F2.cols",header=T,sep="\t")
png_title = "F2_Pigs"; logP_filter = 2; Nsnp = 32197; chrs_sum = 18; y_limit = 17

##################################################################################################
#### [1] Summary the QTLs with -logp>2 of each trait and resort by the trait_class
##################################################################################################

file_all  = dir(); SizeMatrix = length(file_all)
rownames(color_sum) = color_sum[,1]; color_sum_reorder = color_sum[order(color_sum[,2],color_sum[,1]),]
file_all_reorder = color_sum[gsub(".txt.tar.gz","",gsub("Pvalue_","",file_all)),1] # Reorder the phenotypes by the "Trait_Class" and "Trait" #

plotdata_A_all = c(); plotdata_D_all = c(); plotdata_AD_all = c(); plotdata_AvsAD_all = c()
for(i in 1:length(file_all_reorder)){
try({

	tmp_traitname = file_all_reorder[i]
	tmp_filename = paste("Pvalue_",tmp_traitname,".txt.tar.gz",sep="")
	system(paste("tar zxvf ",tmp_filename,sep=""))
	tmp_data = read.table(file=paste(gsub(".tar.gz","",tmp_filename),sep=""),header=T)
	tmp_data = subset(tmp_data,tmp_data[,1]!=0 | tmp_data[,2]!=0)
	tmp_data = tmp_data[complete.cases(tmp_data),]
	system(paste("rm ",gsub(".tar.gz","",tmp_filename),sep=""))
	
	tmp_data_A = cbind(trait=tmp_traitname, subset(tmp_data,tmp_data[,"X.logP.AddCode."] > logP_filter))
	tmp_data_D = cbind(trait=tmp_traitname, subset(tmp_data,tmp_data[,"X.logP.DomCode."] > logP_filter))
	tmp_data_AD = cbind(trait=tmp_traitname, subset(tmp_data,tmp_data[,"X.logP.AddDomCode_Null."] > logP_filter))
	tmp_data_AvsAD = cbind(trait=tmp_traitname, subset(tmp_data,tmp_data[,"X.logP.AddDomCode_Add."] > logP_filter))
	
	plotdata_A_all = rbind(plotdata_A_all,tmp_data_A)
	plotdata_D_all = rbind(plotdata_D_all,tmp_data_D)
	plotdata_AD_all = rbind(plotdata_AD_all,tmp_data_AD)
	plotdata_AvsAD_all = rbind(plotdata_AvsAD_all,tmp_data_AvsAD)

	print(i)
})
}

save(plotdata_A_all, plotdata_D_all, plotdata_AD_all, plotdata_AvsAD_all, color_sum, file=paste("../",png_title,"_SNP_logP2.RData",sep=""))

#################################################################################################################
#### [2] Plot: Integrative manhattan plot using QTLs(-logP>2) from A model, D model, AD model and AvsAD model
#################################################################################################################

### [2-1] Set the plotGWAS function ###
plotGWAS = function(chrs = chrs, plotcolor = c("darkgreen"), alldat = alldat, peakdat = peakdat, y_limit = "", main = "", cex_points = 1.2, cex_lab = 3.1, cex_axis = 3){

	### [STEP1] get the necessary variables ###
    SNP = alldat[,"SNP"]; CHR = alldat[,"chr"]; chrs = 1:chrs
	CHR[CHR == "X"] = chrs[length(chrs)]; CHR[CHR == 23] = chrs[length(chrs)]
    CHR = as.numeric(as.character(CHR))
    POS = as.numeric(alldat[,"pos"])
    idx = which(alldat[,"chr"] %in% chrs)
    if(max(POS,na.rm=T) > 1000){ POS = POS/1e6 } 
	
    chrLen <- tapply(POS,CHR,max); chrmk <- numeric() 
    for(i in 1:length(chrLen)){ chrmk[i] <- sum(as.numeric(chrLen[1:i])) }

    mkrPos = numeric(length(POS)) 
    for(i in as.numeric(unique(CHR))){
        if(i == 1){
            mkrPos[CHR == i] = POS[CHR == i]
        } else {
            mkrPos[CHR == i] = POS[CHR == i] + chrmk[i-1]
        }
    }
    
	### [STEP2] set the plot color ###
	color <- c("darkred","darkgreen","cyan","red", "darkblue", "brown",
               "black","orange","darkred","darkgreen","cyan",
			   "red", "darkblue", "brown","black","orange","darkred",
			   "darkgreen","cyan","red","blue","brown","black","orange")
			   
	# color <- rep(c("#28292B","#C68D34","#31C79C","#0161AD"),6) # Colors from the Qian Li Jiang Shan Plot #
	
    if(length(plotcolor) %in% c(2,3)){ color2 <- rep(plotcolor,length(chrs)) } else { color2 <- color }
	
	### [STEP3] set the plot ylimt ###
	if(y_limit == "NA"){
		y_limit = range(2, 1.1*max(alldat[idx,4],na.rm = TRUE))
	}else{
		y_limit = range(2, 1.1*y_limit)
	}
    
	### [STEP4] start to plot ###
	positions = mkrPos; chridx = which(CHR %in% chrs[1])
	plot(positions[chridx], alldat[chridx, 4],lwd = 1, col = color2[1], 
		 type = "p", pch = 19,ylab = expression(-log[10](italic(P)~value)),xlab ="Chromosomes",
		 xlim = c(0,max(positions,na.rm=T)),main = NULL, xaxt = "n",cex.lab = cex_lab,cex.axis = cex_axis,
		 ylim = y_limit, xaxs = "i", yaxs = "i", cex = cex_points)             
	for(chrom in 2:length(chrs)){
		 chridx = which(CHR %in% chrs[chrom])
		 points(positions[chridx], alldat[chridx, 4], lwd = 1, col = color2[chrom], pch = 19, cex = cex_points)                
	}
	axis(side=1,at=(chrmk-(chrLen/2)),labels=unique(CHR),cex.axis = cex_axis)     

	### [STEP5] add the significant peaks ###
	POS_peak = as.numeric(peakdat[,"pos"]); CHR_peak = peakdat[,"chr"]
	if(max(POS_peak,na.rm=T) > 1000){ POS_peak = POS_peak/1e6 } 
	
	mkrPos_peak = numeric(length(POS_peak)) 
    for(i in as.numeric(unique(CHR_peak))){
        if(i == 1){
            mkrPos_peak[CHR_peak == i] = POS_peak[CHR_peak == i]
        } else {
            mkrPos_peak[CHR_peak == i] = POS_peak[CHR_peak == i] + chrmk[i-1]
        }
    }
	
	positions_peak = mkrPos_peak; chridx_peak = which(CHR_peak %in% chrs[1])
	for(chrom_peak in 1:length(chrs)){
		 chridx_peak = which(CHR_peak %in% chrs[chrom_peak])
		 points(positions_peak[chridx_peak], peakdat[chridx_peak, 4], lwd = 1, col = peakdat[chridx_peak, 5], pch = 18, cex = cex_points*2)                
	}

	return(mkrPos)

}

### [2-2] Get the logP variables of 4 models ###
logP_AddCode = data.frame(SNP=plotdata_A_all[,1], chr=plotdata_A_all[,2], pos=plotdata_A_all[,3], Pc1df1=plotdata_A_all[,4])
logP_DomCode = data.frame(SNP=plotdata_D_all[,1], chr=plotdata_D_all[,2], pos=plotdata_D_all[,3], Pc1df1=plotdata_D_all[,5])
logP_AddDomCode_Null = data.frame(SNP=plotdata_AD_all[,1], chr=plotdata_AD_all[,2], pos=plotdata_AD_all[,3], Pc1df1=plotdata_AD_all[,6])
logP_AddDomCode_Add = data.frame(SNP=plotdata_AvsAD_all[,1], chr=plotdata_AvsAD_all[,2], pos=plotdata_AvsAD_all[,3], Pc1df1=plotdata_AvsAD_all[,7])		

### [2-3] Get the peak significant logP variables of 4 models ###
thrgenome = -log10(0.05/Nsnp); thrsuggest = -log10(1/Nsnp); trait_all = unique(logP_AddCode[,1])

for(i in 1:4){
	all_SigPeak = c()
	logP_all = list(logP_AddCode, logP_DomCode, logP_AddDomCode_Null, logP_AddDomCode_Add)[[i]]
	
	for(j in 1:length(trait_all)){
		# tmp_Sig = subset(logP_all, logP_all[,1]==trait_all[j] & logP_all[,4]>thrgenome)
		tmp_Sig = subset(logP_all, logP_all[,1]==trait_all[j] & logP_all[,4]>thrsuggest)
		
		tmp_SigPeak = c()
		if(nrow(tmp_Sig)>0){
			for(k in 1:length(unique(tmp_Sig[,"chr"]))){
				tmp_Sig_chr = subset(tmp_Sig,tmp_Sig[,"chr"] == unique(tmp_Sig[,"chr"])[k])
				tmp_Sig_max = subset(tmp_Sig_chr,tmp_Sig_chr[,4] == max(tmp_Sig_chr[,4]))[1,]
				tmp_SigPeak = rbind(tmp_SigPeak, tmp_Sig_max)
			}
			all_SigPeak = rbind(all_SigPeak, tmp_SigPeak)
		}
	}
	
	assign(c("logP_AddCode_SigPeak", "logP_DomCode_SigPeak", "logP_AddDomCode_Null_SigPeak", "logP_AddDomCode_Add_SigPeak")[i], all_SigPeak)
}

logP_AddCode_SigPeak = cbind(logP_AddCode_SigPeak, "Color"=color_sum[logP_AddCode_SigPeak[,1],"Color"])
logP_DomCode_SigPeak = cbind(logP_DomCode_SigPeak, "Color"=color_sum[logP_DomCode_SigPeak[,1],"Color"])
logP_AddDomCode_Null_SigPeak = cbind(logP_AddDomCode_Null_SigPeak, "Color"=color_sum[logP_AddDomCode_Null_SigPeak[,1],"Color"])
logP_AddDomCode_Add_SigPeak = cbind(logP_AddDomCode_Add_SigPeak, "Color"=color_sum[logP_AddDomCode_Add_SigPeak[,1],"Color"])

logP_AddCode_SigPeak[logP_AddCode_SigPeak[,4] > y_limit, 4] = y_limit
logP_DomCode_SigPeak[logP_DomCode_SigPeak[,4] > y_limit, 4] = y_limit
logP_AddDomCode_Null_SigPeak[logP_AddDomCode_Null_SigPeak[,4] > y_limit, 4] = y_limit
logP_AddDomCode_Add_SigPeak[logP_AddDomCode_Add_SigPeak[,4] > y_limit, 4] = y_limit

### [2-4] Get the reorder logP variables of 4 models ###
logP_AddCode = logP_AddCode[order(logP_AddCode[,2],logP_AddCode[,3]),]
logP_DomCode = logP_DomCode[order(logP_DomCode[,2],logP_DomCode[,3]),]
logP_AddDomCode_Null = logP_AddDomCode_Null[order(logP_AddDomCode_Null[,2],logP_AddDomCode_Null[,3]),]
logP_AddDomCode_Add = logP_AddDomCode_Add[order(logP_AddDomCode_Add[,2],logP_AddDomCode_Add[,3]),]

logP_AddCode[logP_AddCode[,4] > y_limit, 4] = y_limit
logP_DomCode[logP_DomCode[,4] > y_limit, 4] = y_limit
logP_AddDomCode_Null[logP_AddDomCode_Null[,4] > y_limit, 4] = y_limit
logP_AddDomCode_Add[logP_AddDomCode_Add[,4] > y_limit, 4] = y_limit

### [2-5] Start to plot (horizonal setting or vertical setting) ###
png(paste("../InteManhattan_",png_title,".png",sep=""),width=30000,height=12000,res=600,type="cairo")
layout(matrix(1:4,ncol=2),widths=rep(1,2),heights=rep(1,2))
par(mai=c(1.2,1.5,1,0.5),mgp=c(3.8,1.6,0),cex.main=3.5)

# png(paste("../InteManhattan_",png_title,".png",sep=""),width=15000,height=24000,res=600,type="cairo")
# layout(matrix(1:4,ncol=1),widths=4,heights=rep(1,4))
# par(mai=c(1.2,1.5,1,0.5),mgp=c(3.8,1.6,0),cex.main=3.5)

for(mh_sort in 1:4){
try({
	
	mkrPos = plotGWAS(chrs=chrs_sum, plotcolor=c("#CDB79E","#8B7D6B"), y_limit=y_limit, alldat=na.omit(list(logP_AddCode, logP_AddDomCode_Null, logP_DomCode, logP_AddDomCode_Add)[[mh_sort]]), peakdat=list(logP_AddCode_SigPeak, logP_AddDomCode_Null_SigPeak, logP_DomCode_SigPeak, logP_AddDomCode_Add_SigPeak)[[mh_sort]])
	
	title(paste0(png_title, " : Superimposed Manhattan Plot using ", c("A Model", "AD Model", "D Model", "AvsAD Model")[mh_sort]))
	abline(h=thrgenome,lty=1); abline(h=thrsuggest,lty=2)
	
	if(mh_sort ==4){dev.off()}
})
}
