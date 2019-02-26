### SOURCE ###

#Change source if needed
source1 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/5.1_Subset_diffs_HEGs.csv"
dataH = read.csv(source1, header = TRUE)

source2 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/5.1_Subset_diffs_LEGs.csv"
dataL = read.csv(source2, header = TRUE)

#Column load - LEGs
LAcc <- unlist(dataL[1])
Ltaa.0 <- unlist(dataL[2])
Lnontaa.0 <- unlist(dataL[3])
Ltag.0 <- unlist(dataL[4])
Ltga.0 <- unlist(dataL[5])
Ltaa.1 <- unlist(dataL[6])
Lnontaa.1 <- unlist(dataL[7])
Ltag.1 <- unlist(dataL[8])
Ltga.1 <- unlist(dataL[9])
Ltaa.2 <- unlist(dataL[10])
Lnontaa.2 <- unlist(dataL[11])
Ltag.2 <- unlist(dataL[12])
Ltga.2 <- unlist(dataL[13])
Ltaa.3 <- unlist(dataL[14])
Lnontaa.3 <- unlist(dataL[15])
Ltag.3 <- unlist(dataL[16])
Ltga.3 <- unlist(dataL[17])
Ltaa.4 <- unlist(dataL[18])
Lnontaa.4 <- unlist(dataL[19])
Ltag.4 <- unlist(dataL[20])
Ltga.4 <- unlist(dataL[21])
Ltaa.5 <- unlist(dataL[22])
Lnontaa.5 <- unlist(dataL[23])
Ltag.5 <- unlist(dataL[24])
Ltga.5 <- unlist(dataL[25])
Ltaa.6 <- unlist(dataL[26])
Lnontaa.6 <- unlist(dataL[27])
Ltag.6 <- unlist(dataL[28])
Ltga.6 <- unlist(dataL[29])
Lgc <- unlist(dataL[30])
Lgc3 <- unlist(dataL[31])

#Column load - HEGs
HAcc <- unlist(dataL[1])
Htaa.0 <- unlist(dataL[2])
Hnontaa.0 <- unlist(dataL[3])
Htag.0 <- unlist(dataL[4])
Htga.0 <- unlist(dataL[5])
Htaa.1 <- unlist(dataL[6])
Hnontaa.1 <- unlist(dataL[7])
Htag.1 <- unlist(dataL[8])
Htga.1 <- unlist(dataL[9])
Htaa.2 <- unlist(dataL[10])
Hnontaa.2 <- unlist(dataL[11])
Htag.2 <- unlist(dataL[12])
Htga.2 <- unlist(dataL[13])
Htaa.3 <- unlist(dataL[14])
Hnontaa.3 <- unlist(dataL[15])
Htag.3 <- unlist(dataL[16])
Htga.3 <- unlist(dataL[17])
Htaa.4 <- unlist(dataL[18])
Hnontaa.4 <- unlist(dataL[19])
Htag.4 <- unlist(dataL[20])
Htga.4 <- unlist(dataL[21])
Htaa.5 <- unlist(dataL[22])
Hnontaa.5 <- unlist(dataL[23])
Htag.5 <- unlist(dataL[24])
Htga.5 <- unlist(dataL[25])
Htaa.6 <- unlist(dataL[26])
Hnontaa.6 <- unlist(dataL[27])
Htag.6 <- unlist(dataL[28])
Htga.6 <- unlist(dataL[29])
Hgc <- unlist(dataL[30])
Hgc3 <- unlist(dataL[31])


#Creating dataframes...
total <- merge(dataL, dataH, by="Accession")

pos1 <- total[c("P1_Primary_TAA.x", "P1_Primary_TGA.y")]
rownames(pos1) <- dataL$Accession

pos2 <- total[c("P2_Primary_TAA.x", "P2_Primary_TGA.y")]
rownames(pos2) <- dataL$Accession

pos3 <- total[c("P3_Primary_TAA.x", "P3_Primary_TGA.y")]
rownames(pos3) <- dataL$Accession

pos4 <- total[c("P4_Primary_TAA.x", "P4_Primary_TGA.y")]
rownames(pos4) <- dataL$Accession

pos5 <- total[c("P5_Primary_TAA.x", "P5_Primary_TGA.y")]
rownames(pos5) <- dataL$Accession

pos6 <- total[c("P6_Primary_TAA.x", "P6_Primary_TGA.y")]
rownames(pos6) <- dataL$Accession

#Creating dataframes to contain the differences...
diff1 = (total[["P1_Primary_TGA.y"]] - total[["P1_Primary_TAA.x"]]) / mean(total[["P1_Primary_TGA.y"]] + total[["P1_Primary_TAA.x"]])
diff2 = (total[["P2_Primary_TGA.y"]] - total[["P2_Primary_TAA.x"]]) / mean(total[["P2_Primary_TGA.y"]] + total[["P2_Primary_TAA.x"]])
diff3 = (total[["P3_Primary_TGA.y"]] - total[["P3_Primary_TAA.x"]]) / mean(total[["P3_Primary_TGA.y"]] + total[["P3_Primary_TAA.x"]])
diff4 = (total[["P4_Primary_TGA.y"]] - total[["P4_Primary_TAA.x"]]) / mean(total[["P4_Primary_TGA.y"]] + total[["P4_Primary_TAA.x"]])
diff5 = (total[["P5_Primary_TGA.y"]] - total[["P5_Primary_TAA.x"]]) / mean(total[["P5_Primary_TGA.y"]] + total[["P5_Primary_TAA.x"]])
diff6 = (total[["P6_Primary_TGA.y"]] - total[["P6_Primary_TAA.x"]]) / mean(total[["P6_Primary_TGA.y"]] + total[["P6_Primary_TAA.x"]])

dpos1 = data.frame(diff1, Lgc)
rownames(dpos1) <- dataL$Accession
ordered1 <- dpos1[order(Lgc),] 

dpos2 = data.frame(diff2, Lgc)
rownames(dpos2) <- dataL$Accession
ordered2 <- dpos2[order(Lgc),] 

dpos3 = data.frame(diff3, Lgc)
rownames(dpos3) <- dataL$Accession
ordered3 <- dpos3[order(Lgc),] 

dpos4 = data.frame(diff4, Lgc)
rownames(dpos4) <- dataL$Accession
ordered4 <- dpos4[order(Lgc),] 

dpos5 = data.frame(diff5, Lgc)
rownames(dpos5) <- dataL$Accession
ordered5 <- dpos5[order(Lgc),] 

dpos6 = data.frame(diff6, Lgc)
rownames(dpos6) <- dataL$Accession
ordered6 <- dpos6[order(Lgc),] 

#Plotting differences...

par(mfrow=c(2,3))
barplot(ordered1[["diff1"]], ylab="Standardised frequency difference", ylim=c(-1.5, 2), main="Position +1", cex.main=2, cex.lab=1.5, cex.axis = 1.5, col='#29B6F6')
barplot(ordered2[["diff2"]], ylab="Standardised frequency difference", ylim=c(-1.5, 2), main="Position +2", cex.main=2, cex.lab=1.5, cex.axis = 1.5, col='#29B6F6')
barplot(ordered3[["diff3"]], ylab="Standardised frequency difference", ylim=c(-1.5, 2), main="Position +3", cex.main=2, cex.lab=1.5, cex.axis = 1.5, col='#29B6F6')
barplot(ordered4[["diff4"]], ylab="Standardised frequency difference", ylim=c(-1.5, 2), main="Position +4", cex.main=2, cex.lab=1.5, cex.axis = 1.5, col='#29B6F6')
barplot(ordered5[["diff5"]], ylab="Standardised frequency difference", ylim=c(-1.5, 2), main="Position +5", cex.main=2, cex.lab=1.5, cex.axis = 1.5, col='#29B6F6')
barplot(ordered6[["diff6"]], ylab="Standardised frequency difference", ylim=c(-1.5, 2), main="Position +6", cex.main=2, cex.lab=1.5, cex.axis = 1.5, col='#29B6F6')

#Stats
shapiro.test(Lgc)
shapiro.test(Hgc)
shapiro.test(Lgc3)
shapiro.test(Hgc3)

w1 = wilcox.test(pos1[["P1_Primary_TAA.x"]], pos1[["P1_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
w2 = wilcox.test(pos2[["P2_Primary_TAA.x"]], pos2[["P2_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
w3 = wilcox.test(pos3[["P3_Primary_TAA.x"]], pos3[["P3_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
w4 = wilcox.test(pos4[["P4_Primary_TAA.x"]], pos4[["P4_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
w5 = wilcox.test(pos5[["P5_Primary_TAA.x"]], pos5[["P5_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
w6 = wilcox.test(pos6[["P6_Primary_TAA.x"]], pos6[["P6_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")

library(BSDA)
s1 = SIGN.test(pos1[["P1_Primary_TAA.x"]], pos1[["P1_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
s2 = SIGN.test(pos2[["P2_Primary_TAA.x"]], pos2[["P2_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
s3 = SIGN.test(pos3[["P3_Primary_TAA.x"]], pos3[["P3_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
s4 = SIGN.test(pos4[["P4_Primary_TAA.x"]], pos4[["P4_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
s5 = SIGN.test(pos5[["P5_Primary_TAA.x"]], pos5[["P5_Primary_TGA.y"]], paired = TRUE, alternative = "two.sided")
s6 = SIGN.test(pos6[["P6_Primary_TAA.x"]], pos6[["P6_Primary_TGA.y"]], paired = TRUE, alternative= "two.sided")
