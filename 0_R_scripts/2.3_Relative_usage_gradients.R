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

#Plotting relative codon usage

par(mfrow=c(2,3), cex=1)

plot(gc3, r.taa.1, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.1, col="#fc8d62", pch=16)
points(gc3, r.tga.1, pch=16, col="#8da0cb")
abline(lm(r.taa.1 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.1 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.1 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position 1")

plot(gc3, r.taa.2, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.2, col="#fc8d62", pch=16)
points(gc3, r.tga.2, pch=16, col="#8da0cb")
abline(lm(r.taa.2 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.2 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.2 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position 2")

plot(gc3, r.taa.3, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.3, col="#fc8d62", pch=16)
points(gc3, r.tga.3, pch=16, col="#8da0cb")
abline(lm(r.taa.3 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.3 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.3 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position 3")

plot(gc3, r.taa.4, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.4, col="#fc8d62", pch=16)
points(gc3, r.tga.4, pch=16, col="#8da0cb")
abline(lm(r.taa.4 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.4 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.4 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position 4")

plot(gc3, r.taa.5, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tag.5, col="#fc8d62", pch=16)
points(gc3, r.tga.5, pch=16, col="#8da0cb")
abline(lm(r.taa.5 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.5 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.5 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position 5")

plot(gc3, r.taa.6, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r.tga.6, pch=16, col="#8da0cb")
points(gc3, r.tag.6, col="#fc8d62", pch=16)
abline(lm(r.taa.6 ~ gc3), col="#66c2a5", lwd = 1.5)
abline(lm(r.tag.6 ~ gc3), col="#fc8d62", lwd = 1.5)
abline(lm(r.tga.6 ~ gc3), col="#8da0cb", lwd = 1.5)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("Position 6")

#Get gradients 

lm.taa.0 <- lm(r.taa.0 ~ gc3)
lm.tga.0 <- lm(r.tga.0 ~ gc3)
lm.tag.0 <- lm(r.tag.0 ~ gc3)

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

#Plot gradients
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
sum.1 <- summary(lm.tga.1)
sum.6 <- summary(lm.tga.6)

coeffs.1 <- sum.1$coefficients
coeffs.6 <- sum.6$coefficients

slope.1 <- coeffs.1[2,1]
slope.6 <- coeffs.6[2,1]

sem.1 <- coeffs.1[2,2]
sem.6 <- coeffs.6[2,2]

z <- abs((slope.1-slope.6)/(sqrt(sem.1^2 + sem.6^2)))
pval <- 2*pnorm(z, lower.tail=FALSE)

pval

#Normality testing 
norm.gc <- shapiro.test(gc)     #W = 0.96059, p-value = 3.489e-12
norm.gc3 <- shapiro.test(gc3)   #W = 0.9607, p-value = 3.656e-12

#Correlation testing
spr.stop1.gc3 <- cor.test( ~ stop.1 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9301965 
spr.stop2.gc3 <- cor.test( ~ stop.2 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9494817 
spr.stop3.gc3 <- cor.test( ~ stop.3 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9476725 
spr.stop4.gc3 <- cor.test( ~ stop.4 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9423447 
spr.stop5.gc3 <- cor.test( ~ stop.5 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9391604 
spr.stop6.gc3 <- cor.test( ~ stop.6 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9406134 

spr.taa1.gc3 <- cor.test( ~ r.taa.1 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9396519 
spr.taa2.gc3 <- cor.test( ~ r.taa.2 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9529802 
spr.taa3.gc3 <- cor.test( ~ r.taa.3 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9535181 
spr.taa4.gc3 <- cor.test( ~ r.taa.4 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9475732 
spr.taa5.gc3 <- cor.test( ~ r.taa.5 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9472859 
spr.taa6.gc3 <- cor.test( ~ r.taa.6 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9418209 

spr.tga1.gc3 <- cor.test( ~ r.tga.1 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.6022505 
spr.tga2.gc3 <- cor.test( ~ r.tga.2 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.5828612 
spr.tga3.gc3 <- cor.test( ~ r.tga.3 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.5554547 
spr.tga4.gc3 <- cor.test( ~ r.tga.4 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.4736768 
spr.tga5.gc3 <- cor.test( ~ r.tga.5 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.4518783 
spr.tga6.gc3 <- cor.test( ~ r.tga.6 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.4712456 

spr.tag1.gc3 <- cor.test( ~ r.tag.1 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.7790413 
spr.tag2.gc3 <- cor.test( ~ r.tag.2 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.8226647 
spr.tag3.gc3 <- cor.test( ~ r.tag.3 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.7980667 
spr.tag4.gc3 <- cor.test( ~ r.tag.4 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.7883582 
spr.tag5.gc3 <- cor.test( ~ r.tag.5 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.7693277
spr.tag6.gc3 <- cor.test( ~ r.tag.6 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.7723037 

spr.taa0.gc3 <- cor.test( ~ r.taa.0 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9396519 
spr.tga0.gc3 <- cor.test( ~ r.tga.0 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9396519 
spr.tag0.gc3 <- cor.test( ~ r.tag.0 + gc3, method="spearman") #p-value < 2.2e-16, rho = -0.9396519 

#Obtain best fit equations for TGA at each site
model1 <- lm(tga.1 ~ gc3) 
model2 <- lm(tga.2 ~ gc3)
model3 <- lm(tga.3 ~ gc3) 
model4 <- lm(tga.4 ~ gc3) 
model5 <- lm(tga.5 ~ gc3) 
model6 <- lm(tga.6 ~ gc3) 
