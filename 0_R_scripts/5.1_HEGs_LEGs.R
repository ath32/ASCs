### SOURCE ###

source1 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.1_Frequencies_HEGs.csv"
dataH = read.csv(source1, header = TRUE)

source2 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.1_Frequencies_LEGs.csv"
dataL = read.csv(source2, header = TRUE)

#Loading columns... LEGs
LAcc <- unlist(dataL[1])
Lstop.0 <- unlist(dataL[2])
Ltaa.0 <- unlist(dataL[3])
Ltag.0 <- unlist(dataL[4])
Ltga.0 <- unlist(dataL[5])
Lstop.1 <- unlist(dataL[6])
Ltaa.1 <- unlist(dataL[7])
Ltag.1 <- unlist(dataL[8])
Ltga.1 <- unlist(dataL[9])
Lstop.2 <- unlist(dataL[10])
Ltaa.2 <- unlist(dataL[11])
Ltag.2 <- unlist(dataL[12])
Ltga.2 <- unlist(dataL[13])
Lstop.3 <- unlist(dataL[14])
Ltaa.3 <- unlist(dataL[15])
Ltag.3 <- unlist(dataL[16])
Ltga.3 <- unlist(dataL[17])
Lstop.4 <- unlist(dataL[18])
Ltaa.4 <- unlist(dataL[19])
Ltag.4 <- unlist(dataL[20])
Ltga.4 <- unlist(dataL[21])
Lstop.5 <- unlist(dataL[22])
Ltaa.5 <- unlist(dataL[23])
Ltag.5 <- unlist(dataL[24])
Ltga.5 <- unlist(dataL[25])
Lstop.6 <- unlist(dataL[26])
Ltaa.6<- unlist(dataL[27])
Ltag.6 <- unlist(dataL[28])
Ltga.6 <- unlist(dataL[29])
Lgc <- unlist(dataL[30])
Lgc3 <- unlist(dataL[31])

#Loading columns... HEGs
HAcc <- unlist(dataH[1])
Hstop.0 <- unlist(dataH[2])
Htaa.0 <- unlist(dataH[3])
Htag.0 <- unlist(dataH[4])
Htga.0 <- unlist(dataH[5])
Hstop.1 <- unlist(dataH[6])
Htaa.1 <- unlist(dataH[7])
Htag.1 <- unlist(dataH[8])
Htga.1 <- unlist(dataH[9])
Hstop.2 <- unlist(dataH[10])
Htaa.2 <- unlist(dataH[11])
Htag.2 <- unlist(dataH[12])
Htga.2 <- unlist(dataH[13])
Hstop.3 <- unlist(dataH[14])
Htaa.3 <- unlist(dataH[15])
Htag.3 <- unlist(dataH[16])
Htga.3 <- unlist(dataH[17])
Hstop.4 <- unlist(dataH[18])
Htaa.4 <- unlist(dataH[19])
Htag.4 <- unlist(dataH[20])
Htga.4 <- unlist(dataH[21])
Hstop.5 <- unlist(dataH[22])
Htaa.5 <- unlist(dataH[23])
Htag.5 <- unlist(dataH[24])
Htga.5 <- unlist(dataH[25])
Hstop.6 <- unlist(dataH[26])
Htaa.6<- unlist(dataH[27])
Htag.6 <- unlist(dataH[28])
Htga.6 <- unlist(dataH[29])
Hgc <- unlist(dataH[30])
Hgc3 <- unlist(dataH[31])

#Creating a dataframe that contains all the relevant data
total <- merge(dataL ,dataH, by="Accession")

pos1 <- total[c("P1_stop_f.x", "P1_stop_f.y")]
rownames(pos1) <- dataL$Accession

pos2 <- total[c("P2_stop_f.x", "P2_stop_f.y")]
rownames(pos2) <- dataL$Accession

pos3 <- total[c("P3_stop_f.x", "P3_stop_f.y")]
rownames(pos3) <- dataL$Accession

pos4 <- total[c("P4_stop_f.x", "P4_stop_f.y")]
rownames(pos4) <- dataL$Accession

pos5 <- total[c("P5_stop_f.x", "P5_stop_f.y")]
rownames(pos5) <- dataL$Accession

pos6 <- total[c("P6_stop_f.x", "P6_stop_f.y")]
rownames(pos6) <- dataL$Accession

#Creating dataframes which contain the differences between the two groups
diff1 = total[["P1_stop_f.y"]] - total[["P1_stop_f.x"]]
diff2 = total[["P2_stop_f.y"]] - total[["P1_stop_f.x"]]
diff3 = total[["P3_stop_f.y"]] - total[["P1_stop_f.x"]]
diff4 = total[["P4_stop_f.y"]] - total[["P1_stop_f.x"]]
diff5 = total[["P5_stop_f.y"]] - total[["P1_stop_f.x"]]
diff6 = total[["P6_stop_f.y"]] - total[["P1_stop_f.x"]]

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

#Plotting the differences...

par(mfrow=c(2,3))
barplot(ordered1[["diff1"]], ylab="ASC frequency difference (HEGs - LEGs)", ylim=c(-0.2, 0.1), main="Position +1", cex.main=2, border=NA, cex.lab=1.5, col='#29B6F6')
barplot(ordered2[["diff2"]], ylab="ASC frequency difference (HEGs - LEGs)", ylim=c(-0.2, 0.1), main="Position +2", cex.main=2, border=NA, cex.lab=1.5, col='#29B6F6')
barplot(ordered3[["diff3"]], ylab="ASC frequency difference (HEGs - LEGs)", ylim=c(-0.2, 0.1), main="Position +3", cex.main=2, border=NA, cex.lab=1.5, col='#29B6F6')
barplot(ordered4[["diff4"]], ylab="ASC frequency difference (HEGs - LEGs)", ylim=c(-0.2, 0.1), main="Position +4", cex.main=2, border=NA, cex.lab=1.5, col='#29B6F6')
barplot(ordered5[["diff5"]], ylab="ASC frequency difference (HEGs - LEGs)", ylim=c(-0.2, 0.1), main="Position +5", cex.main=2, border=NA, cex.lab=1.5, col='#29B6F6')
barplot(ordered6[["diff6"]], ylab="ASC frequency difference (HEGs - LEGs)", ylim=c(-0.2, 0.1), main="Position +6", cex.main=2, border=NA, cex.lab=1.5, col='#29B6F6')

#Stats
shapiro.test(Lgc)
shapiro.test(Hgc)
shapiro.test(Lgc3)
shapiro.test(Hgc3)

w1 = wilcox.test(pos1[["P1_stop_f.x"]], pos1[["P1_stop_f.y"]], paired = TRUE, alternative = "two.sided")
w2 = wilcox.test(pos2[["P2_stop_f.x"]], pos2[["P2_stop_f.y"]], paired = TRUE, alternative = "two.sided")
w3 = wilcox.test(pos3[["P3_stop_f.x"]], pos3[["P3_stop_f.y"]], paired = TRUE, alternative = "two.sided")
w4 = wilcox.test(pos4[["P4_stop_f.x"]], pos4[["P4_stop_f.y"]], paired = TRUE, alternative = "two.sided")
w5 = wilcox.test(pos5[["P5_stop_f.x"]], pos5[["P5_stop_f.y"]], paired = TRUE, alternative = "two.sided")
w6 = wilcox.test(pos6[["P6_stop_f.x"]], pos6[["P6_stop_f.y"]], paired = TRUE, alternative = "two.sided")

library(BSDA)
s1 = SIGN.test(pos1[["P1_stop_f.x"]], pos1[["P1_stop_f.y"]], paired = TRUE, alternative = "two.sided")
s2 = SIGN.test(pos2[["P2_stop_f.x"]], pos2[["P2_stop_f.y"]], paired = TRUE, alternative = "two.sided")
s3 = SIGN.test(pos3[["P3_stop_f.x"]], pos3[["P3_stop_f.y"]], paired = TRUE, alternative = "two.sided")
s4 = SIGN.test(pos4[["P4_stop_f.x"]], pos4[["P4_stop_f.y"]], paired = TRUE, alternative = "two.sided")
s5 = SIGN.test(pos5[["P5_stop_f.x"]], pos5[["P5_stop_f.y"]], paired = TRUE, alternative = "two.sided")
s6 = SIGN.test(pos6[["P6_stop_f.x"]], pos6[["P6_stop_f.y"]], paired = TRUE, alternative = "two.sided")


