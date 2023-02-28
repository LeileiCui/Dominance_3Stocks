
rm(list=ls());options(stringsAsFactors=FALSE)

tradir = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE/HSrats_amygdala.phe")
gendir = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes/HSrats_amygdala.phe")
plinkdir = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes/HSrats_amygdala")
plotdat = c("/SAN/mottlab/heterosis/3_HSrats/7_sQTL/Plot_Amygdala/Plot_Amygdala.txt")
bimdat = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes/HSrats_amygdala.bim")

# tradir = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE/HSrats_heart.phe")
# gendir = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes/HSrats_heart.phe")
# plinkdir = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes/HSrats_heart")
# plotdat = c("/SAN/mottlab/heterosis/3_HSrats/7_sQTL/Plot_Heart/Plot_Heart.txt")
# bimdat = c("/SAN/mottlab/heterosis/3_HSrats/1data/RAW_DGE_Genes/HSrats_heart.bim")

library(ggpubr)
dat_tra = read.table(tradir, header=T)
dat_gen = read.table(gendir, header=T)
dat_plot = read.table(plotdat, header=T)
dat_bim = read.table(bimdat, header=F)

for(i in 1:nrow(dat_plot)){

  tmptrait_all = dat_plot[i, c("T1_id","T2_id","G_id","T1_chr","T1_pos")]

  for(j in 1:2){

    try({

      tmptrait = c(as.character(tmptrait_all[c(j, 3)]), subset(dat_bim, dat_bim[,1]==as.numeric(tmptrait_all["T1_chr"]) & dat_bim[,4]==as.numeric(tmptrait_all["T1_pos"]))[,2])
      plot_name = paste(tmptrait[2], tmptrait[1], tmptrait[3], sep="_")

      # [1] Get the expression value of one specified transcript and one specified gene #
      rownames(dat_tra) = paste0("rat", dat_tra[,1])
      rownames(dat_gen) = paste0("rat", dat_gen[,1])
      common_id = intersect(rownames(dat_tra), rownames(dat_gen))
      dat_tra1 = dat_tra[common_id, ]; dat_gen1 = dat_gen[common_id, ]

      # [2] Get the allele value of one specified peak SNP #
      system(paste0("plink --noweb --silent --bfile ", plinkdir, " --snps ", tmptrait[3], " --recode --out tmpsnp"))
      Genoinfo = read.table("tmpsnp.ped", header=F); Genoinfo = Genoinfo[,c(2,7,8)]; system("rm tmpsnp*")
      Genoinfo = cbind("id"=Genoinfo[,1], "allele"=paste(Genoinfo[,2], Genoinfo[,3], sep=""))
      rownames(Genoinfo) = paste0("rat", Genoinfo[,1])
      Genoinfo1 = Genoinfo[common_id, ]

      # [3] Get the final plot data: id, transcript expression, gene expression, allele #
      plotdata = as.data.frame(cbind("id"=1:nrow(dat_tra1), dat_tra1[,tmptrait[1]], dat_gen1[,tmptrait[2]], Genoinfo1[,"allele"]))
      colnames(plotdata) = c("id", tmptrait[1:2], "allele")
      plotdata[,1] = as.numeric(plotdata[,1]); plotdata[,2] = as.numeric(plotdata[,2]); plotdata[,3] = as.numeric(plotdata[,3])

      # [4] Implement the allele QC #
      plotdata = subset(plotdata, plotdata[,"allele"]!="00") # Remove the rows with missing genotype #
      data_clean = plotdata[complete.cases(plotdata),] # Remove the rows with missing value #
      data_clean[,"allele"] = as.factor(data_clean[,"allele"]) # Factorization the allele group #

      # [5] Grouped Scatter plot with marginal density plots #
      png(paste0(plot_name, "_1.png"), width=5000, height=5000, res=600, type="cairo")
      scatter_plot1 <- ggscatterhist(
        data_clean, x = colnames(data_clean)[2], y = colnames(data_clean)[3],
        color = colnames(data_clean)[4], size = 3, alpha = 0.6,
        palette = c("#00AFBB", "#E7B800", "#FC4E07"),
        margin.params = list(fill = colnames(data_clean)[4], color = "black", size = 0.2)
        )
      print(scatter_plot1)
      dev.off()

      # [6] Use box plot as marginal plots #
      png(paste0(plot_name, "_2.png"), width=5000, height=5000, res=600, type="cairo")
      scatter_plot2 <- ggscatterhist(
        data_clean, x = colnames(data_clean)[2], y = colnames(data_clean)[3],
        color = colnames(data_clean)[4], size = 3, alpha = 0.6,
        palette = c("#00AFBB", "#E7B800", "#FC4E07"),
        margin.plot = "boxplot",
        ggtheme = theme_bw()
       )
      print(scatter_plot2)
      dev.off()

    })

    print(paste(i, j, sep=" "))

  }

}

