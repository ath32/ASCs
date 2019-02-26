### SOURCE ### - Choose source file of interest

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/12.1_additional_1_v2.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.6_TAA.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.6_TGA.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.6_TAG.csv"
data = read.csv(source, header = TRUE)

#Calculate differences
diffs <- data[['F_p1']] - data[['F_utr']]
total <- data.frame(diffs, data[['GC3']])
row.names(total) <- data$Accession
colnames(total) <- c('Difference', 'GC3')

#Plot
par(mfrow=c(1,1), cex=1.7)

plot(total[['GC3']], total[['Difference']], xlab='GC3 content (%)', ylab='Raw Frequency Difference', col="#29B6F6",  pch=16, cex.lab=1.5)
abline(h=0, col="black", lwd=3, lty=2)

#Stats

w1 = wilcox.test(data[['F_p1']], data[['F_utr']], paired = TRUE, alternative = "g")

library(BSDA)
s1 = SIGN.test(data[['F_p1']], data[['F_utr']], paired = TRUE, alternative = "g")
