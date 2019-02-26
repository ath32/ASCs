### SOURCE - unhash source file of interest

#Swap source for each plot
source_a <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_allgenes.csv"
source_h <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_HEGs.csv"
source_l <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_LEGs.csv"
#source_tga <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_TGA.csv"
#source_taa <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_TAA.csv"
#source_tag <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_TAG.csv"

data_a = read.csv(source_a, header = TRUE)
data_h = read.csv(source_h, header = TRUE)
data_l = read.csv(source_l, header = TRUE)
#data_taa = read.csv(source_taa, header = TRUE)
#data_tga = read.csv(source_tga, header = TRUE)
#data_tag = read.csv(source_tag, header = TRUE)

#Experimental plots - swap zs for pvalues for different plots

par(mfrow=c(2,3), cex=1.5)

#zscores
plot(data_a[['Genomic_GC3']], data_a[['P1_Zscore']], pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-4,2))
abline(lm.1 <- lm(data_a[['P1_Zscore']]~ data_a[['Genomic_GC3']]), col="black", lwd = 3)
title("All genes (Pos +1)")

plot(data_h[['Genomic_GC3']], data_h[['P1_Zscore']], pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-4,2))
abline(lm.2 <- lm(data_h[['P1_Zscore']]~ data_h[['Genomic_GC3']]), col="black", lwd = 3)
title("HEGs (Pos +1)")

plot(data_l[['Genomic_GC3']], data_l[['P1_Zscore']], pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-4,2))
abline(lm.3 <- lm(data_l[['P1_Zscore']]~ data_l[['Genomic_GC3']]), col="black", lwd = 3)
title("LEGs (Pos +1)")

plot(data_a[['Genomic_GC3']], data_a[['P2_Zscore']], pch=16, col="#01579B", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-4,2))
abline(lm.1 <- lm(data_a[['P2_Zscore']]~ data_a[['Genomic_GC3']]), col="black", lwd = 3)
title("All genes (Pos +2)")

plot(data_h[['Genomic_GC3']], data_h[['P2_Zscore']], pch=16, col="#01579B", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-4,2))
abline(lm.2 <- lm(data_h[['P2_Zscore']]~ data_h[['Genomic_GC3']]), col="black", lwd = 3)
title("HEGs (Pos +2)")

plot(data_l[['Genomic_GC3']], data_l[['P2_Zscore']], pch=16, col="#01579B", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-4,2))
abline(lm.3 <- lm(data_l[['P2_Zscore']]~ data_l[['Genomic_GC3']]), col="black", lwd = 3)
title("LEGs (Pos +2)")

#plot(data_taa[['Genomic_GC3']], data_taa[['P1_Zscore']], pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-5,2))
#title("TAA-terminating genes (Pos +1)")
#abline(lm.3 <- lm(data_taa[['P1_Zscore']]~ data_taa[['Genomic_GC3']]), col="black", lwd = 1.5)

#plot(data_tga[['Genomic_GC3']], data_tga[['P1_Zscore']], pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-5,2))
#abline(lm.4 <- lm(data_tga[['P1_Zscore']]~ data_tga[['Genomic_GC3']]), col="black", lwd = 1.5)
#title("TGA-terminating genes (Pos +1)")

#plot(data_tag[['Genomic_GC3']], data_tag[['P1_Zscore']], pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim = c(-5,2))
#abline(lm.5 <- lm(data_tag[['P1_Zscore']]~ data_tag[['Genomic_GC3']]), col="black", lwd = 1.5)
#title("TAG-terminating genes (Pos +1)")


#Stats for Z-scores
spr.z1.gc3 <- cor.test( ~ data_a[['P1_Zscore']] + data_a[['Genomic_GC3']], method="spearman")
spr.z2.gc3 <- cor.test( ~ data_h[['P1_Zscore']] + data_h[['Genomic_GC3']], method="spearman")
spr.z3.gc3 <- cor.test( ~ data_l[['P1_Zscore']] + data_l[['Genomic_GC3']], method="spearman")
spr.z4.gc3 <- cor.test( ~ data_a[['P2_Zscore']] + data_a[['Genomic_GC3']], method="spearman")
spr.z5.gc3 <- cor.test( ~ data_h[['P2_Zscore']] + data_h[['Genomic_GC3']], method="spearman")
spr.z6.gc3 <- cor.test( ~ data_l[['P2_Zscore']] + data_l[['Genomic_GC3']], method="spearman")
