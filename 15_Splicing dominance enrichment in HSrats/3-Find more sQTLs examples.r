
rm(list=ls());options(stringsAsFactors=FALSE)

Fat_2trpts = read.table("Amygdala_2trpts.txt", header=T)
# Fat_2trpts = read.table("Heart_2trpts.txt", header=T)

Fat_res1_1 = subset(Fat_2trpts, is.na(Fat_2trpts[,"G_logP_AD"]) & !is.na(Fat_2trpts[,"T1_logP_AD"]) & !is.na(Fat_2trpts[,"T2_logP_AD"]) & Fat_2trpts[,"T1_chr"]==Fat_2trpts[,"T2_chr"]) # Find examples of antagonistic T-eQTL, which located in the same chromosome #
Fat_res1_2 = subset(Fat_2trpts, is.na(Fat_2trpts[,"G_logP_AD"]) & !is.na(Fat_2trpts[,"T1_logP_AD"]) & !is.na(Fat_2trpts[,"T2_logP_AD"]) & Fat_2trpts[,"T1_chr"]!=Fat_2trpts[,"T2_chr"]) # Find exmapes of antagonistic T-eQTL, which located in different chromosomes #

Fat_res2_1 = subset(Fat_2trpts, !is.na(Fat_2trpts[,"G_logP_AD"]) & !is.na(Fat_2trpts[,"T1_logP_AD"]) & !is.na(Fat_2trpts[,"T2_logP_AD"]) & Fat_2trpts[,"T1_chr"]==Fat_2trpts[,"T2_chr"] & Fat_2trpts[,"G_chr"]==Fat_2trpts[,"T1_chr"]) # Find examples of synergistic T-eQTL, which located in the same chromosome, and also for the G-eQTL #
Fat_res2_2 = subset(Fat_2trpts, !is.na(Fat_2trpts[,"G_logP_AD"]) & !is.na(Fat_2trpts[,"T1_logP_AD"]) & !is.na(Fat_2trpts[,"T2_logP_AD"]) & Fat_2trpts[,"T1_chr"]==Fat_2trpts[,"T2_chr"] & Fat_2trpts[,"G_chr"]!=Fat_2trpts[,"T1_chr"]) # Find examples of synergistic T-eQTL, which located in the same chromosome, but not for the G-eQTL #
Fat_res2_3 = subset(Fat_2trpts, !is.na(Fat_2trpts[,"G_logP_AD"]) & !is.na(Fat_2trpts[,"T1_logP_AD"]) & !is.na(Fat_2trpts[,"T2_logP_AD"]) & Fat_2trpts[,"T1_chr"]!=Fat_2trpts[,"T2_chr"]) # Find examples of synergistic T-eQTL, which located in different chromosomes, and could be the same as the G-eQTL #

Fat_res = rbind(Fat_res1_1[order(Fat_res1_1[,"T1_logP_AD"],decreasing=T),], 
				colnames(Fat_2trpts), 
				Fat_res1_2[order(Fat_res1_2[,"T1_logP_AD"],decreasing=T),], 
				colnames(Fat_2trpts), 
				Fat_res2_1[order(Fat_res2_1[,"G_logP_AD"],decreasing=T),], 
				colnames(Fat_2trpts), 
				Fat_res2_2[order(Fat_res2_2[,"G_logP_AD"],decreasing=T),], 
				colnames(Fat_2trpts), 
				Fat_res2_3[order(Fat_res2_3[,"G_logP_AD"],decreasing=T),])

write.table(Fat_res, file="Amygdala_examples.txt", row.names=F, col.name=T, quote=F)
# write.table(Fat_res, file="Heart_examples.txt", row.names=F, col.name=T, quote=F)
