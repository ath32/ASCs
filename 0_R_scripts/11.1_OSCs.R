### SOURCE

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/12.1_OSCs.csv"
data = read.csv(source, header = TRUE)

#Calculate relative usage

r0.taa <- data[['X0_TAA']] * 100
r0.tga <- data[['X0_TGA']] * 100
r0.tag <- data[['X0_TAG']] * 100

r1.taa <- data[['X1_TAA']] * 100
r1.tga <- data[['X1_TGA']] * 100
r1.tag <- data[['X1_TAG']] * 100

r2.taa <- data[['X2_TAA']] * 100
r2.tga <- data[['X2_TGA']] * 100
r2.tag <- data[['X2_TAG']] * 100

gc3 <- data[['GC3']]

#Plot 

par(mfrow=c(1,3), cex=1, cex.main=2, cex.lab=1.5, cex.axis=1.5)

plot(gc3, r0.taa, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r0.tag, col="#fc8d62", pch=16)
points(gc3, r0.tga, pch=16, col="#8da0cb")
abline(lm(r0.taa ~ gc3), col="#66c2a5", lwd = 3)
abline(lm(r0.tag ~ gc3), col="#fc8d62", lwd = 3)
abline(lm(r0.tga ~ gc3), col="#8da0cb", lwd = 3)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("In-frame")

plot(gc3, r1.taa, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r1.tag, col="#fc8d62", pch=16)
points(gc3, r1.tga, pch=16, col="#8da0cb")
abline(lm(r1.taa ~ gc3), col="#66c2a5", lwd = 3)
abline(lm(r1.tag ~ gc3), col="#fc8d62", lwd = 3)
abline(lm(r1.tga ~ gc3), col="#8da0cb", lwd = 3)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("+1 Frame-shift")

plot(gc3, r2.taa, pch=16, col="#66c2a5", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)")
points(gc3, r2.tag, col="#fc8d62", pch=16)
points(gc3, r2.tga, pch=16, col="#8da0cb")
abline(lm(r2.taa ~ gc3), col="#66c2a5", lwd = 3)
abline(lm(r2.tag ~ gc3), col="#fc8d62", lwd = 3)
abline(lm(r2.tga ~ gc3), col="#8da0cb", lwd = 3)
legend("topright", pch=16, col = c("#66c2a5", "#fc8d62", "#8da0cb"), legend = c("TAA", "TAG", "TGA"))
title("+2 Frame-shift")


#Spearman's ranks

s0_taa <- cor.test( ~ r0.taa + gc3, method="spearman")
s0_tga <- cor.test( ~ r0.tga + gc3, method="spearman")
s0_tag <- cor.test( ~ r0.tag + gc3, method="spearman")

s1_taa <- cor.test( ~ r1.taa + gc3, method="spearman")
s1_tga <- cor.test( ~ r1.tga + gc3, method="spearman")
s1_tag <- cor.test( ~ r1.tag + gc3, method="spearman")

s2_taa <- cor.test( ~ r2.taa + gc3, method="spearman")
s2_tga <- cor.test( ~ r2.tga + gc3, method="spearman")
s2_tag <- cor.test( ~ r2.tag + gc3, method="spearman")




