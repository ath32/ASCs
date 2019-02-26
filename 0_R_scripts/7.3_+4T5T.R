### SOURCE ###

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.2_+4T5T.csv"
data = read.csv(source, header = TRUE)

#Calculate differences
diffs <- data[['F_4T5T_genes']] - data[['F_TT_codons']]
total <- data.frame(diffs, data[['GC3']])
row.names(total) <- data$Accession
colnames(total) <- c('Difference', 'GC3')

#Plot
par(mfrow=c(1,1), cex=1.5)

plot(total[['GC3']], total[['Difference']], xlab='GC3 Content (%)', ylab='Difference', col="#29B6F6",  pch=16)
abline(h=0, col="black", lwd=3, lty=2)

#Stats

w1 = wilcox.test(data[['F_4T5T_genes']], data[['F_TT_codons']], paired = TRUE, alternative = "two.sided")

library(BSDA)
s1 = SIGN.test(data[['F_4T5T_genes']], data[['F_TT_codons']], paired = TRUE, alternative = "two.sided")
