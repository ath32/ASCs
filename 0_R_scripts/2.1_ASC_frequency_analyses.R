### SOURCE ###

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.1_Frequencies_all.csv"
data = read.csv(source, header = TRUE)

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

#Calculate relative codon usage

r.taa.0 <- taa.0 / stop.0 * 100
r.tga.0 <- tga.0 / stop.0 * 100
r.tag.0 <- tag.0 / stop.0 * 100

r.taa.1 <- taa.1 / stop.1 * 100
r.tga.1 <- tga.1 / stop.1 * 100
r.tag.1 <- tag.1 / stop.1 * 100

r.taa.2 <- taa.2 / stop.2 * 100
r.tga.2 <- tga.2 / stop.2 * 100
r.tag.2 <- tag.2 / stop.2 * 100

r.taa.3 <- taa.3 / stop.3 * 100
r.tga.3 <- tga.3 / stop.3 * 100
r.tag.3 <- tag.3 / stop.3 * 100

r.taa.4 <- taa.4 / stop.4 * 100
r.tga.4 <- tga.4 / stop.4 * 100
r.tag.4 <- tag.4 / stop.4 * 100

r.taa.5 <- taa.5 / stop.5 * 100
r.tga.5 <- tga.5 / stop.5 * 100
r.tag.5 <- tag.5 / stop.5 * 100

r.taa.6 <- taa.6 / stop.6 * 100
r.tga.6 <- tga.6 / stop.6 * 100
r.tag.6 <- tag.6 / stop.6 * 100

r.taa.0 <- taa.0 / stop.0 * 100
r.tga.0 <- tga.0 / stop.0 * 100
r.tag.0 <- tag.0 / stop.0 * 100

#Plot 1 - Stops at each position vs gc3

par(mfrow=c(2,2), cex.lab=2, cex.axis=1.5)

plot(gc3, stop.1, pch=16, cex=0.7, col="#1b9e77", xlab = "GC3 content (%)", ylab = "Stop codon frequency")
points(gc3, stop.2, col="#d95f02", pch=16, cex=0.7)
points(gc3, stop.3, pch=16, col="#7570b3", cex=0.7)
points(gc3, stop.4, pch=16, col="#e7298a", cex=0.7)
points(gc3, stop.5, pch=16, col="#66a61e", cex=0.7)
points(gc3, stop.6, pch=16, col="#e6ab02", cex=0.7)
abline(g1.stop.1 <- lm(stop.1 ~ gc3), col="#1b9e77", lwd = 1.5)
abline(g1.stop.2 <- lm(stop.2 ~ gc3), col="#d95f02", lwd = 1.5)
abline(g1.stop.3 <- lm(stop.3 ~ gc3), col="#7570b3", lwd = 1.5)
abline(g1.stop.4 <- lm(stop.4 ~ gc3), col="#e7298a", lwd = 1.5)
abline(g1.stop.5 <- lm(stop.5 ~ gc3), col="#66a61e", lwd = 1.5)
abline(g1.stop.6 <- lm(stop.6 ~ gc3), col="#e6ab02", lwd = 1.5)
legend("topright", pch=16, col = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"), cex=1.2, legend = c("position 1", "position 2", "position 3", "position 4", "position 5", "position 6"))

plot(gc3, taa.1, pch=16, cex=0.7, col="#1b9e77", xlab = "GC3 content (%)", ylab = "TAA codon frequency")
points(gc3, taa.2, col="#d95f02", pch=16, cex=0.7)
points(gc3, taa.3, pch=16, col="#7570b3", cex=0.7)
points(gc3, taa.4, pch=16, col="#e7298a", cex=0.7)
points(gc3, taa.5, pch=16, col="#66a61e", cex=0.7)
points(gc3, taa.6, pch=16, col="#e6ab02", cex=0.7)
abline(g1.taa.1 <- lm(taa.1 ~ gc3), col="#1b9e77", lwd = 1.5)
abline(g1.taa.2 <- lm(taa.2 ~ gc3), col="#d95f02", lwd = 1.5)
abline(g1.taa.3 <- lm(taa.3 ~ gc3), col="#e7298a", lwd = 1.5)
abline(g1.taa.4 <- lm(taa.4 ~ gc3), col="#7570b3", lwd = 1.5)
abline(g1.taa.5 <- lm(taa.5 ~ gc3), col="#66a61e", lwd = 1.5)
abline(g1.taa.6 <- lm(taa.6 ~ gc3), col="#e6ab02", lwd = 1.5)
legend("topright", pch=16, col = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"), cex=1.2, legend = c("position 1", "position 2", "position 3", "position 4", "position 5", "position 6"))

plot(gc3, tag.1, pch=16, cex=0.7, col="#1b9e77", xlab = "GC3 content (%)", ylab = "TAG codon frequency")
points(gc3, tag.2, col="#d95f02", pch=16, cex=0.7)
points(gc3, tag.3, pch=16, col="#7570b3", cex=0.7)
points(gc3, tag.4, pch=16, col="#e7298a", cex=0.7)
points(gc3, tag.5, pch=16, col="#66a61e", cex=0.7)
points(gc3, tag.6, pch=16, col="#e6ab02", cex=0.7)
abline(g1.tag.1 <- lm(tag.1 ~ gc3), col="#1b9e77", lwd = 1.5)
abline(g1.tag.2 <- lm(tag.2 ~ gc3), col="#7570b3", lwd = 1.5)
abline(g1.tag.3 <- lm(tag.3 ~ gc3), col="#d95f02", lwd = 1.5)
abline(g1.tag.4 <- lm(tag.4 ~ gc3), col="#e7298a", lwd = 1.5)
abline(g1.tag.5 <- lm(tag.5 ~ gc3), col="#66a61e", lwd = 1.5)
abline(g1.tag.6 <- lm(tag.6 ~ gc3), col="#e6ab02", lwd = 1.5)
legend("topright", pch=16, col = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"), cex=1.2, legend = c("position 1", "position 2", "position 3", "position 4", "position 5", "position 6"))

plot(gc3, tga.1, pch=16, cex=0.7, col="#1b9e77", xlab = "GC3 content (%)", ylab = "TGA codon frequency")
points(gc3, tga.2, col="#d95f02", pch=16, cex=0.7)
points(gc3, tga.3, pch=16, col="#7570b3", cex=0.7)
points(gc3, tga.4, pch=16, col="#e7298a", cex=0.7)
points(gc3, tga.5, pch=16, col="#66a61e", cex=0.7)
points(gc3, tga.6, pch=16, col="#e6ab02", cex=0.7)
abline(g1.tga.1 <- lm(tga.1 ~ gc3), col="#1b9e77", lwd = 1.5)
abline(g1.tga.2 <- lm(tga.2 ~ gc3), col="#d95f02", lwd = 1.5)
abline(g1.tga.3 <- lm(tga.3 ~ gc3), col="#7570b3", lwd = 1.5)
abline(g1.tga.4 <- lm(tga.4 ~ gc3), col="#e7298a", lwd = 1.5)
abline(g1.tga.5 <- lm(tga.5 ~ gc3), col="#66a61e", lwd = 1.5)
abline(g1.tga.6 <- lm(tga.6 ~ gc3), col="#e6ab02", lwd = 1.5)
legend("topright", pch=16, col = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"), cex=1.2, legend = c("position 1", "position 2", "position 3", "position 4", "position 5", "position 6"))

#Get gradients for above plot
coef(g1.stop.1)[2] 
coef(g1.stop.2)[2] 
coef(g1.stop.3)[2] 
coef(g1.stop.4)[2] 
coef(g1.stop.5)[2] 
coef(g1.stop.6)[2] 

coef(g1.taa.1)[2] 
coef(g1.taa.2)[2] 
coef(g1.taa.3)[2] 
coef(g1.taa.4)[2] 
coef(g1.taa.5)[2] 
coef(g1.taa.6)[2] 

coef(g1.tag.1)[2] 
coef(g1.tag.2)[2] 
coef(g1.tag.3)[2] 
coef(g1.tag.4)[2] 
coef(g1.tag.5)[2] 
coef(g1.tag.6)[2] 

coef(g1.tga.1)[2] 
coef(g1.tga.2)[2] 
coef(g1.tga.3)[2] 
coef(g1.tga.4)[2] 
coef(g1.tga.5)[2] 
coef(g1.tga.6)[2] 

#Plot gradients
par(mfrow=c(2,2), cex=1)

stop_bar_1 <- barplot(c(coef(g1.stop.1)[2], coef(g1.stop.2)[2], coef(g1.stop.3)[2], coef(g1.stop.4)[2], coef(g1.stop.5)[2], coef(g1.stop.6)[2]), ylim=c(-0.003,0.002), names.arg=c("1", "2", "3", "4", "5", "6"), ylab = "Gradient", xlab = "Position", main="Any stop")
lines(x=stop_bar_1, y=c(coef(g1.stop.1)[2], coef(g1.stop.2)[2], coef(g1.stop.3)[2], coef(g1.stop.4)[2], coef(g1.stop.5)[2], coef(g1.stop.6)[2]))
points(x=stop_bar_1, y=c(coef(g1.stop.1)[2], coef(g1.stop.2)[2], coef(g1.stop.3)[2], coef(g1.stop.4)[2], coef(g1.stop.5)[2], coef(g1.stop.6)[2]))

taa_bar_1 <- barplot(c(coef(g1.taa.1)[2], coef(g1.taa.2)[2], coef(g1.taa.3)[2], coef(g1.taa.4)[2], coef(g1.taa.5)[2], coef(g1.taa.6)[2]), ylim=c(-0.003,0.002), names.arg=c("1", "2", "3", "4", "5", "6"), ylab = "Gradient", xlab = "Position", main="TAA")
lines(x=taa_bar_1, y=c(coef(g1.taa.1)[2], coef(g1.taa.2)[2], coef(g1.taa.3)[2], coef(g1.taa.4)[2], coef(g1.taa.5)[2], coef(g1.taa.6)[2]))
points(x=taa_bar_1, y=c(coef(g1.taa.1)[2], coef(g1.taa.2)[2], coef(g1.taa.3)[2], coef(g1.taa.4)[2], coef(g1.taa.5)[2], coef(g1.taa.6)[2]))

tga_bar_1 <- barplot(c(coef(g1.tga.1)[2], coef(g1.tga.2)[2], coef(g1.tga.3)[2], coef(g1.tga.4)[2], coef(g1.tga.5)[2], coef(g1.tga.6)[2]), ylim=c(-0.0005,0.0005), names.arg=c("1", "2", "3", "4", "5", "6"), ylab = "Gradient", xlab = "Position", main="TGA")
lines(x=tga_bar_1, y=c(coef(g1.tga.1)[2], coef(g1.tga.2)[2], coef(g1.tga.3)[2], coef(g1.tga.4)[2], coef(g1.tga.5)[2], coef(g1.tga.6)[2]))
points(x=tga_bar_1, y=c(coef(g1.tga.1)[2], coef(g1.tga.2)[2], coef(g1.tga.3)[2], coef(g1.tga.4)[2], coef(g1.tga.5)[2], coef(g1.tga.6)[2]))

tag_bar_1 <- barplot(c(coef(g1.tag.1)[2], coef(g1.tag.2)[2], coef(g1.tag.3)[2], coef(g1.tag.4)[2], coef(g1.tag.5)[2], coef(g1.tag.6)[2]), ylim=c(-0.0005,0.0005), names.arg=c("1", "2", "3", "4", "5", "6"), ylab = "Gradient", xlab = "Position", main="TAG")
lines(x=tag_bar_1, y=c(coef(g1.tag.1)[2], coef(g1.tag.2)[2], coef(g1.tag.3)[2], coef(g1.tag.4)[2], coef(g1.tag.5)[2], coef(g1.tag.6)[2]))
points(x=tag_bar_1, y=c(coef(g1.tag.1)[2], coef(g1.tag.2)[2], coef(g1.tag.3)[2], coef(g1.tag.4)[2], coef(g1.tag.5)[2], coef(g1.tag.6)[2]))

#Are gradients at pos1 and pos6 different?
sum.1 <- summary(g1.tag.1)
sum.6 <- summary(g1.tag.6)

coeffs.1 <- sum.1$coefficients
coeffs.6 <- sum.6$coefficients

slope.1 <- coeffs.1[2,1]
slope.6 <- coeffs.6[2,1]

sem.1 <- coeffs.1[2,2]
sem.6 <- coeffs.6[2,2]

z <- abs((slope.1-slope.6)/(sqrt(sem.1^2 + sem.6^2)))
pval <- 2*pnorm(z, lower.tail=FALSE)

#Plotting mean frequencies at each position 

par(mfrow=c(2,2), cex=0.7)

barplot(c(mean(stop.1), mean(stop.2), mean(stop.3), mean(stop.4), mean(stop.5), mean(stop.6)), width = 1, beside = FALSE, ylim=c(0,0.06), xlab="Codon Position", ylab='Observed Frequency', names.arg=c("1", "2", "3", "4", "5", "6"), main = "Frequency of stop codons at each 3' position")
abline(h=0.046875)

barplot(c(mean(taa.1), mean(taa.2), mean(taa.3), mean(taa.4), mean(taa.5), mean(taa.6)), width = 1, beside = FALSE, ylim=c(0,0.03), xlab="Codon Position", ylab='Observed Frequency', names.arg=c("1", "2", "3", "4", "5", "6"), main = "Frequency of taa codons at each 3' position")
abline(h=0.015625)

barplot(c(mean(tag.1), mean(tag.2), mean(tag.3), mean(tag.4), mean(tag.5), mean(tag.6)), width = 1, beside = FALSE, ylim=c(0,0.03), xlab="Codon Position", ylab='Observed Frequency', names.arg=c("1", "2", "3", "4", "5", "6"), main = "Frequency of tag codons at each 3' position")
abline(h=0.015625)

barplot(c(mean(tga.1), mean(tga.2), mean(tga.3), mean(tga.4), mean(tga.5), mean(tga.6)), width = 1, beside = FALSE, ylim=c(0,0.03), xlab="Codon Position", ylab='Observed Frequency', names.arg=c("1", "2", "3", "4", "5", "6"), main = "Frequency of tga codons at each 3' position")
abline(h=0.015625)

Plotting stop codon usage at each position

par(mfrow=c(2,3), cex=1)

plot(gc3, r.taa.1, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.1, col="#fc8d62", pch=16)
points(gc3, r.tga.1, pch=16, col="#8da0cb")
abline(lm(r.taa.1 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.1 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.1 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position +1")

plot(gc3, r.taa.2, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.2, col="#fc8d62", pch=16)
points(gc3, r.tga.2, pch=16, col="#8da0cb")
abline(lm(r.taa.2 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.2 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.2 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position +2")

plot(gc3, r.taa.3, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.3, col="#fc8d62", pch=16)
points(gc3, r.tga.3, pch=16, col="#8da0cb")
abline(lm(r.taa.3 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.3 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.3 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position +3")

plot(gc3, r.taa.4, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.4, col="#fc8d62", pch=16)
points(gc3, r.tga.4, pch=16, col="#8da0cb")
abline(lm(r.taa.4 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.4 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.4 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position +4")

plot(gc3, r.taa.5, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.5, col="#fc8d62", pch=16)
points(gc3, r.tga.5, pch=16, col="#8da0cb")
abline(lm(r.taa.5 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.5 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.5 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position +5")

plot(gc3, r.taa.6, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.6, col="#fc8d62", pch=16)
points(gc3, r.tga.6, pch=16, col="#8da0cb")
abline(lm(r.taa.6 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.6 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.6 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position +6")

#Get gradients
lm.taa.1 <- lm(r.taa.1 ~ gc3)
lm.tga.1 <- lm(r.tga.1 ~ gc3)
lm.tag.1 <- lm(r.tag.1 ~ gc3)

lm.taa.2 <- lm(r.taa.2 ~ gc3)
lm.tga.2 <- lm(r.tga.2 ~ gc3)
lm.tag.2 <- lm(r.tag.2 ~ gc3)

lm.taa.3 <- lm(r.taa.3 ~ gc3)
lm.tga.3 <- lm(r.tga.3 ~ gc3)
lm.tag.3 <- lm(r.tag.3 ~ gc3)

lm.taa.4 <- lm(r.taa.4 ~ gc3)
lm.tga.4 <- lm(r.tga.4 ~ gc3)
lm.tag.4 <- lm(r.tag.4 ~ gc3)

lm.taa.5 <- lm(r.taa.5 ~ gc3)
lm.tga.5 <- lm(r.tga.5 ~ gc3)
lm.tag.5 <- lm(r.tag.5 ~ gc3)

lm.taa.6 <- lm(r.taa.6 ~ gc3)
lm.tga.6 <- lm(r.tga.6 ~ gc3)
lm.tag.6 <- lm(r.tag.6 ~ gc3)

g.taa.1 <- coef(lm.taa.1)[2]
g.tga.1 <- coef(lm.tga.1)[2]
g.tag.1 <- coef(lm.tag.1)[2]

g.taa.2 <- coef(lm.taa.2)[2]
g.tga.2 <- coef(lm.tga.2)[2]
g.tag.2 <- coef(lm.tag.2)[2]

g.taa.3 <- coef(lm.taa.3)[2]
g.tga.3 <- coef(lm.tga.3)[2]
g.tag.3 <- coef(lm.tag.3)[2]

g.taa.4 <- coef(lm.taa.4)[2]
g.tga.4 <- coef(lm.tga.4)[2]
g.tag.4 <- coef(lm.tag.4)[2]

g.taa.5 <- coef(lm.taa.5)[2]
g.tga.5 <- coef(lm.tga.5)[2]
g.tag.5 <- coef(lm.tag.5)[2]

g.taa.6 <- coef(lm.taa.6)[2]
g.tga.6 <- coef(lm.tga.6)[2]
g.tag.6 <- coef(lm.tag.6)[2]

#Plotting gradients
par(mfrow=c(1,3), cex=1)

taa_bar <- barplot(c(g.taa.1, g.taa.2, g.taa.3, g.taa.4, g.taa.5, g.taa.6), ylim=c(-1.5,1.5), names.arg=c("1", "2", "3", "4", "5", "6"), ylab = "Gradient", xlab = "Position", main="TAA")
lines(x=taa_bar, y=c(g.taa.1, g.taa.2, g.taa.3, g.taa.4, g.taa.5, g.taa.6))
points(x=taa_bar, y=c(g.taa.1, g.taa.2, g.taa.3, g.taa.4, g.taa.5, g.taa.6))

tga_bar <- barplot(c(g.tga.1, g.tga.2, g.tga.3, g.tga.4, g.tga.5, g.tga.6), ylim=c(-1.5,1.5), names.arg=c("1", "2", "3", "4", "5", "6"), ylab = "Gradient", xlab = "Position", main="TGA")
lines(x=tga_bar, y=c(g.tga.1, g.tga.2, g.tga.3, g.tga.4, g.tga.5, g.tga.6))
points(x=tga_bar, y=c(g.tga.1, g.tga.2, g.tga.3, g.tga.4, g.tga.5, g.tga.6))

tag_bar <- barplot(c(g.tag.1, g.tag.2, g.tag.3, g.tag.4, g.tag.5, g.tag.6), ylim=c(-1.5,1.5), names.arg=c("1", "2", "3", "4", "5", "6"), ylab = "Gradient", xlab = "Position", main="TAG")
lines(x=tag_bar, y=c(g.tag.1, g.tag.2, g.tag.3, g.tag.4, g.tag.5, g.tag.6))
points(x=tag_bar, y=c(g.tag.1, g.tag.2, g.tag.3, g.tag.4, g.tag.5, g.tag.6))

#Are gradients at pos1 and pos6 different?
sum.1 <- summary(lm.tag.1)
sum.6 <- summary(lm.tag.6)

coeffs.1 <- sum.1$coefficients
coeffs.6 <- sum.6$coefficients

slope.1 <- coeffs.1[2,1]
slope.6 <- coeffs.6[2,1]

sem.1 <- coeffs.1[2,2]
sem.6 <- coeffs.6[2,2]

z <- abs((slope.1-slope.6)/(sqrt(sem.1^2 + sem.6^2)))
pval <- 2*pnorm(z, lower.tail=FALSE)
slope.1
slope.6
pval

#Normality testing 
norm.gc <- shapiro.test(gc)     #W = 0.96059, p-value = 3.489e-12
norm.gc3 <- shapiro.test(gc3)   #W = 0.9607, p-value = 3.656e-12

#Correlation testing
spr.stop1.gc3 <- cor.test( ~ stop.1 + gc3, method="spearman") 
spr.stop2.gc3 <- cor.test( ~ stop.2 + gc3, method="spearman") 
spr.stop3.gc3 <- cor.test( ~ stop.3 + gc3, method="spearman") 
spr.stop4.gc3 <- cor.test( ~ stop.4 + gc3, method="spearman") 
spr.stop5.gc3 <- cor.test( ~ stop.5 + gc3, method="spearman") 
spr.stop6.gc3 <- cor.test( ~ stop.6 + gc3, method="spearman") 

spr.taa1.gc3 <- cor.test( ~ taa.1 + gc3, method="spearman") 
spr.taa2.gc3 <- cor.test( ~ taa.2 + gc3, method="spearman")
spr.taa3.gc3 <- cor.test( ~ taa.3 + gc3, method="spearman") 
spr.taa4.gc3 <- cor.test( ~ taa.4 + gc3, method="spearman")
spr.taa5.gc3 <- cor.test( ~ taa.5 + gc3, method="spearman")
spr.taa6.gc3 <- cor.test( ~ taa.6 + gc3, method="spearman") 

spr.tga1.gc3 <- cor.test( ~ tga.1 + gc3, method="spearman") 
spr.tga2.gc3 <- cor.test( ~ tga.2 + gc3, method="spearman") 
spr.tga3.gc3 <- cor.test( ~ tga.3 + gc3, method="spearman") 
spr.tga4.gc3 <- cor.test( ~ tga.4 + gc3, method="spearman") 
spr.tga5.gc3 <- cor.test( ~ tga.5 + gc3, method="spearman")
spr.tga6.gc3 <- cor.test( ~ tga.6 + gc3, method="spearman") 

spr.tag1.gc3 <- cor.test( ~ tag.1 + gc3, method="spearman")
spr.tag2.gc3 <- cor.test( ~ tag.2 + gc3, method="spearman")
spr.tag3.gc3 <- cor.test( ~ tag.3 + gc3, method="spearman")
spr.tag4.gc3 <- cor.test( ~ tag.4 + gc3, method="spearman")
spr.tag5.gc3 <- cor.test( ~ tag.5 + gc3, method="spearman")
spr.tag6.gc3 <- cor.test( ~ tag.6 + gc3, method="spearman")

#Obtain best fit equations for TGA at each site
model1 <- lm(tga.1 ~ gc3) 
model2 <- lm(tga.2 ~ gc3) 
model3 <- lm(tga.3 ~ gc3) 
model4 <- lm(tga.4 ~ gc3) 
model5 <- lm(tga.5 ~ gc3)
model6 <- lm(tga.6 ~ gc3)





