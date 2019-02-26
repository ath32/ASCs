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

#Plotting stop codon usage at each position

par(mfrow=c(2,4), cex=1)

plot(gc3, r.taa.0, pch=16, col="#B3F5FC", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)", ylim=c(0,100))
points(gc3, r.tag.0, col="#01579B", pch=16)
points(gc3, r.tga.0, pch=16, col="#29B6F6")
legend("topright", pch=16, col = c("#B3F5FC", "#01579B", "#29B6F6"), legend = c("TAA", "TAG", "TGA"), cex=1)
title("Primary Stop (Position +0)")

plot(gc3, r.taa.1, pch=16, col="#B3F5FC", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)", ylim=c(0,100))
points(gc3, r.tag.1, col="#01579B", pch=16)
points(gc3, r.tga.1, pch=16, col="#29B6F6")
legend("topright", pch=16, col = c("#B3F5FC", "#01579B", "#29B6F6"), legend = c("TAA", "TAG", "TGA"))
title("Position +1")

plot(gc3, r.taa.2, pch=16, col="#B3F5FC", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)", ylim=c(0,100))
points(gc3, r.tag.2, col="#01579B", pch=16)
points(gc3, r.tga.2, pch=16, col="#29B6F6")
legend("topright", pch=16, col = c("#B3F5FC", "#01579B", "#29B6F6"), legend = c("TAA", "TAG", "TGA"))
title("Position +2")

plot(gc3, r.taa.3, pch=16, col="#B3F5FC", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)", ylim=c(0,100))
points(gc3, r.tag.3, col="#01579B", pch=16)
points(gc3, r.tga.3, pch=16, col="#29B6F6")
legend("topright", pch=16, col = c("#B3F5FC", "#01579B", "#29B6F6"), legend = c("TAA", "TAG", "TGA"))
title("Position +3")

plot(gc3, r.taa.4, pch=16, col="#B3F5FC", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)", ylim=c(0,100))
points(gc3, r.tag.4, col="#01579B", pch=16)
points(gc3, r.tga.4, pch=16, col="#29B6F6")
legend("topright", pch=16, col = c("#B3F5FC", "#01579B", "#29B6F6"), legend = c("TAA", "TAG", "TGA"))
title("Position +4")

plot(gc3, r.taa.5, pch=16, col="#B3F5FC", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)", ylim=c(0,100))
points(gc3, r.tag.5, col="#01579B", pch=16)
points(gc3, r.tga.5, pch=16, col="#29B6F6")
legend("topright", pch=16, col = c("#B3F5FC", "#01579B", "#29B6F6"), legend = c("TAA", "TAG", "TGA"))
title("Position +5")

plot(gc3, r.taa.6, pch=16, col="#B3F5FC", xlab = "GC3 content (%)", ylab = "Stop codon usage (%)", ylim=c(0,100))
points(gc3, r.tag.6, col="#01579B", pch=16)
points(gc3, r.tga.6, pch=16, col="#29B6F6")
legend("topright", pch=16, col = c("#B3F5FC", "#01579B", "#29B6F6"), legend = c("TAA", "TAG", "TGA"))
title("Position +6")

#Correlation testing
shapiro.test(gc)
shapiro.test(gc3)

spr.r.taa.0 <- cor.test( ~ r.taa.0 + gc3, method="spearman") 
spr.r.tag.0 <- cor.test( ~ r.tag.0 + gc3, method="spearman") 
spr.r.tga.0 <- cor.test( ~ r.tga.0 + gc3, method="spearman") 

spr.r.taa.1 <- cor.test( ~ r.taa.1 + gc3, method="spearman") 
spr.r.tag.1 <- cor.test( ~ r.tag.1 + gc3, method="spearman") 
spr.r.tga.1 <- cor.test( ~ r.tga.1 + gc3, method="spearman") 

spr.r.taa.2 <- cor.test( ~ r.taa.2 + gc3, method="spearman") 
spr.r.tag.2 <- cor.test( ~ r.tag.2 + gc3, method="spearman") 
spr.r.tga.2 <- cor.test( ~ r.tga.2 + gc3, method="spearman") 

spr.r.taa.3 <- cor.test( ~ r.taa.3 + gc3, method="spearman") 
spr.r.tag.3 <- cor.test( ~ r.tag.3 + gc3, method="spearman") 
spr.r.tga.3 <- cor.test( ~ r.tga.3 + gc3, method="spearman") 

spr.r.taa.4 <- cor.test( ~ r.taa.4 + gc3, method="spearman") 
spr.r.tag.4 <- cor.test( ~ r.tag.4 + gc3, method="spearman") 
spr.r.tga.4 <- cor.test( ~ r.tga.4 + gc3, method="spearman") 

spr.r.taa.5 <- cor.test( ~ r.taa.5 + gc3, method="spearman") 
spr.r.tag.5 <- cor.test( ~ r.tag.5 + gc3, method="spearman") 
spr.r.tga.5 <- cor.test( ~ r.tga.5 + gc3, method="spearman") 

spr.r.taa.6 <- cor.test( ~ r.taa.6 + gc3, method="spearman") 
spr.r.tag.6 <- cor.test( ~ r.tag.6 + gc3, method="spearman") 
spr.r.tga.6 <- cor.test( ~ r.tga.6 + gc3, method="spearman") 

