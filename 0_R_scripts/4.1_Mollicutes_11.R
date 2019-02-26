### SOURCE ###

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.1_Frequencies_all.csv"
data = read.csv(source, header = TRUE)

observed_4 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/4.2_TT11_mollicute_ASCs.csv"
data_4 = read.csv(observed_4, header = TRUE)

observed_11 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/4.2_TT11_mollicute_ASCs.csv"
data_11 = read.csv(observed_11, header = TRUE)

#Loading columns...
Acc <- unlist(data[1])
stop.0 <- unlist(data[2])
taa.0 <- unlist(data[3])
tag.0 <- unlist(data[4])
tga.0 <- unlist(data[5])
stop.1 <- unlist(data[6])
taa.1 <- unlist(data[7])
tag.1 <- unlist(data[8])
tga.1 <- unlist(data[9])
stop.2 <- unlist(data[10])
taa.2 <- unlist(data[11])
tag.2 <- unlist(data[12])
tga.2 <- unlist(data[13])
stop.3 <- unlist(data[14])
taa.3 <- unlist(data[15])
tag.3 <- unlist(data[16])
tga.3 <- unlist(data[17])
stop.4 <- unlist(data[18])
taa.4 <- unlist(data[19])
tag.4 <- unlist(data[20])
tga.4 <- unlist(data[21])
stop.5 <- unlist(data[22])
taa.5 <- unlist(data[23])
tag.5 <- unlist(data[24])
tga.5 <- unlist(data[25])
stop.6 <- unlist(data[26])
taa.6<- unlist(data[27])
tag.6 <- unlist(data[28])
tga.6 <- unlist(data[29])
gc <- unlist(data[30])
gc3 <- unlist(data[31])

#Plot 1 - Stops at each position vs gc3

par(mfrow=c(2,3), cex=1.3, cex.main=2)

scatter.smooth(gc3, tga.1, pch=16, col="#1b9e77", xlab = "GC3 content (%)", ylab = "TGA codon frequency", main='Position +1', ylim=c(0,0.045))
loess.smooth(gc3, tga.1, span = 2/3, degree = 1, col="#1b9e77")

scatter.smooth(gc3, tga.2, pch=16, col="#d95f02", xlab = "GC3 content (%)", ylab = "TGA codon frequency", main='Position +2', ylim=c(0,0.045))
loess.smooth(gc3, tga.2, span = 2/3, degree = 1, col="#d95f02")

scatter.smooth(gc3, tga.3, pch=16, col="#7570b3", xlab = "GC3 content (%)", ylab = "TGA codon frequency", main='Position +3', ylim=c(0,0.045))
loess.smooth(gc3, tga.3, span = 2/3, degree = 1, col="#d7570b3")

scatter.smooth(gc3, tga.4, pch=16, col="#e7298a", xlab = "GC3 content (%)", ylab = "TGA codon frequency", main='Position +4', ylim=c(0,0.045))
loess.smooth(gc3, tga.4, span = 2/3, degree = 1, col="#e7298a")

scatter.smooth(gc3, tga.5, pch=16, col="#66a61e", xlab = "GC3 content (%)", ylab = "TGA codon frequency", main='Position +5', ylim=c(0,0.045))
loess.smooth(gc3, tga.5, span = 2/3, degree = 1, col="#66a61e")

scatter.smooth(gc3, tga.6, pch=16, col="#e6ab02", xlab = "GC3 content (%)", ylab = "TGA codon frequency", main='Position +6', ylim=c(0,0.045))
loess.smooth(gc3, tga.6, span = 2/3, degree = 1, col="#e6ab02")

#Predictions
p1 <- c()
p2 <- c()
p3 <- c()
p4 <- c()
p5 <- c()
p6 <- c()

for (row in 1:nrow(data_4)) {
  
  o.1 <- data_4[row, "P1_Observed"]
  o.2 <- data_4[row, "P2_Observed"]
  o.3 <- data_4[row, "P3_Observed"]
  o.4 <- data_4[row, "P4_Observed"]
  o.5 <- data_4[row, "P5_Observed"]
  o.6 <- data_4[row, "P6_Observed"]
  
  p.1 <- predict(object = loess(tga.1 ~ gc3, span=2/3, control=loess.control(surface="direct")), newdata = data_4[row, "Genomic_GC3"])
  p.2 <- predict(object = loess(tga.2 ~ gc3, span=2/3, control=loess.control(surface="direct")), newdata = data_4[row, "Genomic_GC3"])
  p.3 <- predict(object = loess(tga.3 ~ gc3, span=2/3, control=loess.control(surface="direct")), newdata = data_4[row, "Genomic_GC3"])
  p.4 <- predict(object = loess(tga.4 ~ gc3, span=2/3, control=loess.control(surface="direct")), newdata = data_4[row, "Genomic_GC3"])
  p.5 <- predict(object = loess(tga.5 ~ gc3, span=2/3), control=loess.control(surface="direct"), newdata = data_4[row, "Genomic_GC3"])
  p.6 <- predict(object = loess(tga.6 ~ gc3, span=2/3, control=loess.control(surface="direct")), newdata = data_4[row, "Genomic_GC3"])
  
  p1 <- c(p1, p.1)
  p2 <- c(p2, p.2)
  p3 <- c(p3, p.3)
  p4 <- c(p4, p.4)
  p5 <- c(p5, p.5)
  p6 <- c(p6, p.6)
  
}


#Create dataframe of new values and append to original dataset
d1 <- data.frame(p1)
rownames(d1) <- data_4[['Accession']]
total4 <- cbind(data_4, d1)

d2 <- data.frame(p2)
rownames(d2) <- data_4[['Accession']]
total4 <- cbind(total4, d2)

d3 <- data.frame(p3)
rownames(d3) <- data_4[['Accession']]
total4 <- cbind(total4, d3)

d4 <- data.frame(p4)
rownames(d4) <- data_4[['Accession']]
total4 <- cbind(total4, d4)

d5 <- data.frame(p5)
rownames(d5) <- data_4[['Accession']]
total4 <- cbind(total4, d5)

d6 <- data.frame(p6)
rownames(d6) <- data_4[['Accession']]
total4 <- cbind(total4, d6)

#Create columns for differences
total4['Diff1'] <- total4['P1_Observed'] - total4['p1']
total4['Diff2'] <- total4['P2_Observed'] - total4['p2']
total4['Diff3'] <- total4['P3_Observed'] - total4['p3']
total4['Diff4'] <- total4['P4_Observed'] - total4['p4']
total4['Diff5'] <- total4['P5_Observed'] - total4['p5']
total4['Diff6'] <- total4['P6_Observed'] - total4['p6']

#Plotting bars of difference between observed and predicted

par(mfrow=c(2,3), cex=1)

barplot(total4[['Diff1']], ylab = 'Observed - Predicted', ylim=c(-0.02, 0.04), main='Pos +1', col='lightblue')
barplot(total4[['Diff2']], ylab = 'Observed - Predicted', ylim=c(-0.02, 0.04), main='Pos +2', col='lightblue')
barplot(total4[['Diff3']], ylab = 'Observed - Predicted', ylim=c(-0.02, 0.04), main='Pos +3', col='lightblue')
barplot(total4[['Diff4']], ylab = 'Observed - Predicted', ylim=c(-0.02, 0.04), main='Pos +4', col='lightblue')
barplot(total4[['Diff5']], ylab = 'Observed - Predicted', ylim=c(-0.02, 0.04), main='Pos +5', col='lightblue')
barplot(total4[['Diff6']], ylab = 'Observed - Predicted', ylim=c(-0.02, 0.04), main='Pos +6', col='lightblue')

#Plotting scatter of difference against GC3 content 

par(mfrow=c(2,3), cex=1)

plot(total4[['GC3']], total4[['Diff1']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
abline(h=0, lty=2, col='red', lwd=3)
title("Position +1")

plot(total4[['GC3']], total4[['Diff2']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
abline(h=0, lty=2, col='red', lwd=3)
title("Position +2")

plot(total4[['GC3']], total4[['Diff3']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
abline(h=0, lty=2, col='red', lwd=3)
title("Position +3")

plot(total4[['GC3']], total4[['Diff4']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
abline(h=0, lty=2, col='red', lwd=3)
title("Position +4")

plot(total4[['GC3']], total4[['Diff5']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
abline(h=0, lty=2, col='red', lwd=3)
title("Position +5")

plot(total4[['GC3']], total4[['Diff6']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
abline(h=0, lty=2, col='red', lwd=3)
title("Position +6")


#Statistics to compare observed to expected
w1 = wilcox.test(total4[['P1_Observed']], total4[['p1']], paired = TRUE, alternative = "less")
w2 = wilcox.test(total4[['P2_Observed']], total4[['p2']], paired = TRUE, alternative = "less")
w3 = wilcox.test(total4[['P3_Observed']], total4[['p3']], paired = TRUE, alternative = "less")
w4 = wilcox.test(total4[['P4_Observed']], total4[['p4']], paired = TRUE, alternative = "less")
w5 = wilcox.test(total4[['P5_Observed']], total4[['p5']], paired = TRUE, alternative = "less")
w6 = wilcox.test(total4[['P6_Observed']], total4[['p6']], paired = TRUE, alternative = "less")

library(BSDA)
s1 = SIGN.test(total4[['P1_Observed']], total4[['p1']], paired = TRUE, alternative = "less")
s2 = SIGN.test(total4[['P2_Observed']], total4[['p2']], paired = TRUE, alternative = "less")
s3 = SIGN.test(total4[['P3_Observed']], total4[['p3']], paired = TRUE, alternative = "less")
s4 = SIGN.test(total4[['P4_Observed']], total4[['p4']], paired = TRUE, alternative = "less")
s5 = SIGN.test(total4[['P5_Observed']], total4[['p5']], paired = TRUE, alternative = "less")
s6 = SIGN.test(total4[['P6_Observed']], total4[['p6']], paired = TRUE, alternative = "less")


