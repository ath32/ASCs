### SOURCE ### 

source_all <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.3_Freq_pri_all.csv"
data_all = read.csv(source_all, header = TRUE)

source_hegs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.3_Freq_pri_HEGs.csv"
data_hegs = read.csv(source_hegs, header = TRUE)

source_legs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.3_Freq_pri_LEGs.csv"
data_legs = read.csv(source_legs, header = TRUE)

#Get means for each codon for each position - all
a.taa.1 <- mean(data_all[['TAA_1']])
a.taa.2 <- mean(data_all[['TAA_2']])
a.taa.3 <- mean(data_all[['TAA_3']])
a.taa.4 <- mean(data_all[['TAA_4']])
a.taa.5 <- mean(data_all[['TAA_5']])
a.taa.6 <- mean(data_all[['TAA_6']])

a.tga.1 <- mean(data_all[['TGA_1']])
a.tga.2 <- mean(data_all[['TGA_2']])
a.tga.3 <- mean(data_all[['TGA_3']])
a.tga.4 <- mean(data_all[['TGA_4']])
a.tga.5 <- mean(data_all[['TGA_5']])
a.tga.6 <- mean(data_all[['TGA_6']])

a.tag.1 <- mean(data_all[['TAG_1']])
a.tag.2 <- mean(data_all[['TAG_2']])
a.tag.3 <- mean(data_all[['TAG_3']])
a.tag.4 <- mean(data_all[['TAG_4']])
a.tag.5 <- mean(data_all[['TAG_5']])
a.tag.6 <- mean(data_all[['TAG_6']])

#Get means for each codon for each position - hegs
h.taa.1 <- mean(data_hegs[['TAA_1']])
h.taa.2 <- mean(data_hegs[['TAA_2']])
h.taa.3 <- mean(data_hegs[['TAA_3']])
h.taa.4 <- mean(data_hegs[['TAA_4']])
h.taa.5 <- mean(data_hegs[['TAA_5']])
h.taa.6 <- mean(data_hegs[['TAA_6']])

h.tga.1 <- mean(data_hegs[['TGA_1']])
h.tga.2 <- mean(data_hegs[['TGA_2']])
h.tga.3 <- mean(data_hegs[['TGA_3']])
h.tga.4 <- mean(data_hegs[['TGA_4']])
h.tga.5 <- mean(data_hegs[['TGA_5']])
h.tga.6 <- mean(data_hegs[['TGA_6']])

h.tag.1 <- mean(data_hegs[['TAG_1']])
h.tag.2 <- mean(data_hegs[['TAG_2']])
h.tag.3 <- mean(data_hegs[['TAG_3']])
h.tag.4 <- mean(data_hegs[['TAG_4']])
h.tag.5 <- mean(data_hegs[['TAG_5']])
h.tag.6 <- mean(data_hegs[['TAG_6']])

#Get means for each codon for each position - legs
l.taa.1 <- mean(data_legs[['TAA_1']])
l.taa.2 <- mean(data_legs[['TAA_2']])
l.taa.3 <- mean(data_legs[['TAA_3']])
l.taa.4 <- mean(data_legs[['TAA_4']])
l.taa.5 <- mean(data_legs[['TAA_5']])
l.taa.6 <- mean(data_legs[['TAA_6']])

l.tga.1 <- mean(data_legs[['TGA_1']])
l.tga.2 <- mean(data_legs[['TGA_2']])
l.tga.3 <- mean(data_legs[['TGA_3']])
l.tga.4 <- mean(data_legs[['TGA_4']])
l.tga.5 <- mean(data_legs[['TGA_5']])
l.tga.6 <- mean(data_legs[['TGA_6']])

l.tag.1 <- mean(data_legs[['TAG_1']])
l.tag.2 <- mean(data_legs[['TAG_2']])
l.tag.3 <- mean(data_legs[['TAG_3']])
l.tag.4 <- mean(data_legs[['TAG_4']])
l.tag.5 <- mean(data_legs[['TAG_5']])
l.tag.6 <- mean(data_legs[['TAG_6']])

#Create dataframes for each position
all.1 <- c(a.taa.1, a.tga.1, a.tag.1)
all.2 <- c(a.taa.2, a.tga.2, a.tag.2)
all.3 <- c(a.taa.3, a.tga.3, a.tag.3)
all.4 <- c(a.taa.4, a.tga.4, a.tag.4)
all.5 <- c(a.taa.5, a.tga.5, a.tag.5)
all.6 <- c(a.taa.6, a.tga.6, a.tag.6)

a.total <- rbind(all.1, all.2, all.3, all.4, all.5, all.6)
colnames(a.total) <- c('TAA', 'TGA', 'TAG')
rownames(a.total) <- c('+1', '+2', '+3', '+4', '+5', '+6')

hegs.1 <- c(h.taa.1, h.tga.1, h.tag.1)
hegs.2 <- c(h.taa.2, h.tga.2, h.tag.2)
hegs.3 <- c(h.taa.3, h.tga.3, h.tag.3)
hegs.4 <- c(h.taa.4, h.tga.4, h.tag.4)
hegs.5 <- c(h.taa.5, h.tga.5, h.tag.5)
hegs.6 <- c(h.taa.6, h.tga.6, h.tag.6)

h.total <- rbind(hegs.1, hegs.2, hegs.3, hegs.4, hegs.5, hegs.6)
colnames(h.total) <- c('TAA', 'TGA', 'TAG')
rownames(h.total) <- c('+1', '+2', '+3', '+4', '+5', '+6')

legs.1 <- c(l.taa.1, l.tga.1, l.tag.1)
legs.2 <- c(l.taa.2, l.tga.2, l.tag.2)
legs.3 <- c(l.taa.3, l.tga.3, l.tag.3)
legs.4 <- c(l.taa.4, l.tga.4, l.tag.4)
legs.5 <- c(l.taa.5, l.tga.5, l.tag.5)
legs.6 <- c(l.taa.6, l.tga.6, l.tag.6)

l.total <- rbind(legs.1, legs.2, legs.3, legs.4, legs.5, legs.6)
colnames(l.total) <- c('TAA', 'TGA', 'TAG')
rownames(l.total) <- c('+1', '+2', '+3', '+4', '+5', '+6')

#Plot barplot
par(mfrow=c(1,3), cex=1.5)

barplot(t(as.matrix((a.total))), main='ALL', beside=TRUE, col=c("#B3F5FC", "#29B6F6", "#01579B"), ylab="Frequency", ylim=c(0,0.12), cex.lab=1.5, border=NA, xlab = 'Position', cex.axis=1.3, cex.names=1.3, cex.main = 1.5)
legend("topright", 
       legend = c("TAA-terminating", "TGA-terminating", "TAG-terminating"),
       fill = c("#B3F5FC", "#29B6F6", "#01579B"), bty = "n")

barplot(t(as.matrix((h.total))), main='HEGs', beside=TRUE, col=c("#B3F5FC", "#29B6F6", "#01579B"), ylab="Frequency", ylim=c(0,0.12), cex.lab=1.5, border=NA, xlab = 'Position', cex.axis=1.3, cex.names=1.3, cex.main = 1.5)
legend("topright", 
       legend = c("TAA-terminating", "TGA-terminating", "TAG-terminating"),
       fill = c("#B3F5FC", "#29B6F6", "#01579B"), bty = "n")

barplot(t(as.matrix((l.total))), main='LEGs', beside=TRUE, col=c("#B3F5FC", "#29B6F6", "#01579B"), ylab="Frequency", ylim=c(0,0.12), cex.lab=1.5, border=NA, xlab = 'Position', cex.axis=1.3, cex.names=1.3, cex.main = 1.5)
legend("topright", 
       legend = c("TAA-terminating", "TGA-terminating", "TAG-terminating"),
       fill = c("#B3F5FC", "#29B6F6", "#01579B"), bty = "n")

#Various statistics
a1.taa = wilcox.test(data_all[['TGA_1']], data_all[['TAA_1']], paired = TRUE, alternative = "greater")
a1.tag = wilcox.test(data_all[['TGA_1']], data_all[['TAG_1']], paired = TRUE, alternative = "greater")

a2.taa = wilcox.test(data_all[['TGA_2']], data_all[['TAA_2']], paired = TRUE, alternative = "greater")
a2.tag = wilcox.test(data_all[['TGA_2']], data_all[['TAG_2']], paired = TRUE, alternative = "greater")

a3.taa = wilcox.test(data_all[['TGA_3']], data_all[['TAA_3']], paired = TRUE, alternative = "greater")
a3.tag = wilcox.test(data_all[['TGA_3']], data_all[['TAG_3']], paired = TRUE, alternative = "greater")

a4.taa = wilcox.test(data_all[['TGA_4']], data_all[['TAA_4']], paired = TRUE, alternative = "greater")
a4.tag = wilcox.test(data_all[['TGA_4']], data_all[['TAG_4']], paired = TRUE, alternative = "greater")

a5.taa = wilcox.test(data_all[['TGA_5']], data_all[['TAA_5']], paired = TRUE, alternative = "greater")
a5.tag = wilcox.test(data_all[['TGA_5']], data_all[['TAG_5']], paired = TRUE, alternative = "greater")

a6.taa = wilcox.test(data_all[['TGA_6']], data_all[['TAA_6']], paired = TRUE, alternative = "greater")
a6.tag = wilcox.test(data_all[['TGA_6']], data_all[['TAG_6']], paired = TRUE, alternative = "greater")

a <- kruskal.test(list(data_legs[['TAA_1']], data_legs[['TGA_1']], data_legs[['TAG_1']]))
b <- kruskal.test(list(data_hegs[['TAA_2']], data_hegs[['TGA_2']], data_hegs[['TAG_2']]))
c <- kruskal.test(list(data_legs[['TAA_3']], data_legs[['TGA_3']], data_legs[['TAG_3']]))
d <- kruskal.test(list(data_legs[['TAA_4']], data_legs[['TGA_4']], data_legs[['TAG_4']]))
e <- kruskal.test(list(data_legs[['TAA_5']], data_legs[['TGA_5']], data_legs[['TAG_5']]))
f <- kruskal.test(list(data_legs[['TAA_6']], data_legs[['TGA_6']], data_legs[['TAG_6']]))
g <- wilcox.test(data_hegs[['TGA_2']], data_hegs[['TAA_2']], paired = TRUE, alternative = "greater")
h <- wilcox.test(data_hegs[['TGA_2']], data_hegs[['TAG_2']], paired = TRUE, alternative = "greater")
j <- wilcox.test(data_hegs[['TGA_1']], data_legs[['TGA_1']], paired = FALSE, alternative = "greater")



