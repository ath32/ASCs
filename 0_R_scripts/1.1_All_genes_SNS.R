### SOURCE ###

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/2.2_Simulations.csv"
data = read.csv(source, header = TRUE)

#Loading columns...
Acc <- unlist(data[1])
ob.0 <- unlist(data[2])
tot.0 <- unlist(data[3])
exf.0 <- unlist(data[4])
obf.0 <- unlist(data[5])
sd.0 <- unlist(data[6])
z.0 <- unlist(data[7])
binomP.0 <- unlist(data[8])
ob.1 <- unlist(data[9])
tot.1 <- unlist(data[10])
exf.1 <- unlist(data[11])
obf.1 <- unlist(data[12])
sd.1 <- unlist(data[13])
z.1 <- unlist(data[14])
binomP.1 <- unlist(data[15])
ob.2 <- unlist(data[16])
tot.2 <- unlist(data[17])
exf.2 <- unlist(data[18])
obf.2 <- unlist(data[19])
sd.2 <- unlist(data[20])
z.2 <- unlist(data[21])
binomP.2 <- unlist(data[22])
ob.3 <- unlist(data[23])
tot.3 <- unlist(data[24])
exf.3 <- unlist(data[25])
obf.3 <- unlist(data[26])
sd.3 <- unlist(data[27])
z.3 <- unlist(data[28])
binomP.3 <- unlist(data[29])
ob.4 <- unlist(data[30])
tot.4 <- unlist(data[31])
exf.4 <- unlist(data[32])
obf.4 <- unlist(data[33])
sd.4 <- unlist(data[34])
z.4 <- unlist(data[35])
binomP.4 <- unlist(data[36])
ob.5 <- unlist(data[37])
tot.5 <- unlist(data[38])
exf.5 <- unlist(data[39])
obf.5 <- unlist(data[40])
sd.5 <- unlist(data[41])
z.5 <- unlist(data[42])
binomP.5 <- unlist(data[43])
ob.6 <- unlist(data[44])
tot.6 <- unlist(data[45])
exf.6 <- unlist(data[46])
obf.6 <- unlist(data[47])
sd.6 <- unlist(data[48])
z.6 <- unlist(data[49])
binomP.6 <- unlist(data[50])
gc <- unlist(data[51])
gc3 <- unlist(data[52])

#Experimental plot 1 - P-values

par(mfrow=c(1,1), cex=1.2)

plot(gc3, binomP.1, pch=16, col="#fc8d62", xlab = "GC3 content (%)", ylab = "Two-tailed Binomial Test P-value")
title("Position 1")

plot(gc3, binomP.2, pch=16, col="#fc8d62", xlab = "GC3 content (%)", ylab = "Two-tailed Binomial Test P-value")
title("Position 2")

plot(gc3, binomP.3, pch=16, col="#fc8d62", xlab = "GC3 content (%)", ylab = "Two-tailed Binomial Test P-value")
title("Position 3")

plot(gc3, binomP.4, pch=16, col="#fc8d62", xlab = "GC3 content (%)", ylab = "Two-tailed Binomial Test P-value")
title("Position 4")

plot(gc3, binomP.5, pch=16, col="#fc8d62", xlab = "GC3 content (%)", ylab = "Two-tailed Binomial Test P-value")
title("Position 5")

plot(gc3, binomP.6, pch=16, col="#fc8d62", xlab = "GC3 content (%)", ylab = "Two-tailed Binomial Test P-value")
title("Position 6")

#Experimental plot 1 - Z-scores

par(mfrow=c(1,1), cex=1.2)

plot(gc3, z.1, pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim=c(-8, 2.5), cex.lab=1.5, cex.axis = 1.5)
abline(lm.z1 <- lm(z.1 ~ gc3), col="black", lwd = 2.5)
title("Position +1", cex.main=2.5)

plot(gc3, z.2, pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim=c(-8, 2.5), cex.lab=2, cex.axis = 1.5)
abline(lm.z2 <- lm(z.2 ~ gc3), col="black", lwd = 2.5)
title("Position +2", cex.main=2.5)

plot(gc3, z.3, pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim=c(-8, 2.5), cex.lab=2, cex.axis = 1.5)
abline(lm.z3 <- lm(z.3 ~ gc3), col="black", lwd = 2.5)
title("Position +3", cex.main=2.5)

plot(gc3, z.4, pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim=c(-8, 2.5), cex.lab=2, cex.axis = 1.5)
abline(lm.z4 <- lm(z.4 ~ gc3), col="black", lwd = 2.5)
title("Position +4", cex.main=2.5)

plot(gc3, z.5, pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim=c(-8, 2.5), cex.lab=2, cex.axis = 1.5)
abline(lm.z5 <- lm(z.5 ~ gc3), col="black", lwd = 2.5)
title("Position +5", cex.main=2.5)

plot(gc3, z.6, pch=16, col="#29B6F6", xlab = "GC3 content (%)", ylab = "Z-score", ylim=c(-8, 2.5), cex.lab=2, cex.axis = 1.5)
abline(lm.z6 <- lm(z.6 ~ gc3), col="black", lwd = 2.5)
title("Position +6", cex.main=2.5)

#Spearman's ranks for Z scores
spr.z1.gc3 <- cor.test( ~ z.1 + gc3, method="spearman") 
spr.z2.gc3 <- cor.test( ~ z.2 + gc3, method="spearman") 
spr.z3.gc3 <- cor.test( ~ z.3 + gc3, method="spearman")
spr.z4.gc3 <- cor.test( ~ z.4 + gc3, method="spearman") 
spr.z5.gc3 <- cor.test( ~ z.5 + gc3, method="spearman") 
spr.z6.gc3 <- cor.test( ~ z.6 + gc3, method="spearman")

#Calculate gradients of lines
g.z1 <- coef(lm.z1)[2]
g.z2 <- coef(lm.z2)[2]
g.z3 <- coef(lm.z3)[2]
g.z4 <- coef(lm.z4)[2]
g.z5 <- coef(lm.z5)[2]
g.z6 <- coef(lm.z6)[2]

#Are gradients at posX and posY different? (sub in different linear models for different comparisons)
sum.x <- summary(lm.z5)
sum.y <- summary(lm.z6)

coeffs.x <- sum.x$coefficients
coeffs.y <- sum.y$coefficients

slope.x <- coeffs.x[2,1]
slope.y <- coeffs.y[2,1]

sem.x <- coeffs.x[2,2]
sem.y <- coeffs.y[2,2]

z <- abs((slope.x-slope.y)/(sqrt(sem.x^2 + sem.y^2)))
pval <- 2*pnorm(z, lower.tail=FALSE)

#Wilcoxons - sub in different columns for different comparisons
w1 = wilcox.test(obf.6, exf.6, paired = TRUE, alternative = "greater")
w2 = wilcox.test(obf.6, exf.6, paired = TRUE, alternative = "l")
