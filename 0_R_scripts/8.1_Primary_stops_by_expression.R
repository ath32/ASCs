### SOURCE

source1 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/8.1_primary_stop_percentages_HEGs.csv"
dataH = read.csv(source1, header = TRUE)

source2 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/8.1_primary_stop_percentages_LEGs.csv"
dataL = read.csv(source2, header = TRUE)

#Calculate means
mean.taa.h <- mean(dataH[['TAA.']])
mean.tga.h <- mean(dataH[['TGA.']])
mean.tag.h <- mean(dataH[['TAG.']])

mean.taa.l <- mean(dataL[['TAA.']])
mean.tga.l <- mean(dataL[['TGA.']])
mean.tag.l <- mean(dataL[['TAG.']])

#Create vectors
hegs <- c(mean.taa.h, mean.tga.h, mean.tag.h)
legs <- c(mean.taa.l, mean.tga.l, mean.tag.l)
labels <- c('TAA', 'TGA', 'TAG')

#Plot pies

par(mfrow=c(1,2), mar = c(5,0,4,0), cex=1.5)
pie(hegs, labels, col=c("#99d594", "#fc8d59", 'yellow'))
title("HEGs", line = -1)
pie(legs, labels, col=c("#99d594", "#fc8d59", 'yellow'))
title("LEGs", line = -1)

#Stats
taa <- wilcox.test(dataH[['TAA.']], dataL[['TAA.']], paired = TRUE, alternative = "greater")
tga <- wilcox.test(dataH[['TGA.']], dataL[['TGA.']], paired = TRUE, alternative = "l")
tag <- wilcox.test(dataH[['TAG.']], dataL[['TAG.']], paired = TRUE, alternative = "l")