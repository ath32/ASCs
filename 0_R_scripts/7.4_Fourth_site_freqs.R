### SOURCE ###

#Change source if needed
source_all <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.4_fourth_site_all.csv"
data_all = read.csv(source_all, header = TRUE)

source_HEGs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.4_fourth_site_hegs.csv"
data_HEGs = read.csv(source_HEGs, header = TRUE)

source_LEGs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.4_fourth_site_legs.csv"
data_LEGs = read.csv(source_LEGs, header = TRUE)

#Calculate means - ALL
taa_a <- mean(data_all[['TAA_FourthA']])
taa_t <- mean(data_all[['TAA_FourthT']])
taa_g <- mean(data_all[['TAA_FourthG']])
taa_c <- mean(data_all[['TAA_FourthC']])
tga_a <- mean(data_all[['TGA_FourthA']])
tga_t <- mean(data_all[['TGA_FourthT']])
tga_g <- mean(data_all[['TGA_FourthG']])
tga_c <- mean(data_all[['TGA_FourthC']])
tag_a <- mean(data_all[['TAG_FourthA']])
tag_t <- mean(data_all[['TAG_FourthT']])
tag_g <- mean(data_all[['TAG_FourthG']])
tag_c <- mean(data_all[['TAG_FourthC']])

#Calculate means - HEGs
HEGs_taa_a <- mean(data_HEGs[['TAA_FourthA']])
HEGs_taa_t <- mean(data_HEGs[['TAA_FourthT']])
HEGs_taa_g <- mean(data_HEGs[['TAA_FourthG']])
HEGs_taa_c <- mean(data_HEGs[['TAA_FourthC']])
HEGs_tga_a <- mean(data_HEGs[['TGA_FourthA']])
HEGs_tga_t <- mean(data_HEGs[['TGA_FourthT']])
HEGs_tga_g <- mean(data_HEGs[['TGA_FourthG']])
HEGs_tga_c <- mean(data_HEGs[['TGA_FourthC']])
HEGs_tag_a <- mean(data_HEGs[['TAG_FourthA']])
HEGs_tag_t <- mean(data_HEGs[['TAG_FourthT']])
HEGs_tag_g <- mean(data_HEGs[['TAG_FourthG']])
HEGs_tag_c <- mean(data_HEGs[['TAG_FourthC']])

#Calculate means - LEGs
LEGs_taa_a <- mean(data_LEGs[['TAA_FourthA']])
LEGs_taa_t <- mean(data_LEGs[['TAA_FourthT']])
LEGs_taa_g <- mean(data_LEGs[['TAA_FourthG']])
LEGs_taa_c <- mean(data_LEGs[['TAA_FourthC']])
LEGs_tga_a <- mean(data_LEGs[['TGA_FourthA']])
LEGs_tga_t <- mean(data_LEGs[['TGA_FourthT']])
LEGs_tga_g <- mean(data_LEGs[['TGA_FourthG']])
LEGs_tga_c <- mean(data_LEGs[['TGA_FourthC']])
LEGs_tag_a <- mean(data_LEGs[['TAG_FourthA']])
LEGs_tag_t <- mean(data_LEGs[['TAG_FourthT']])
LEGs_tag_g <- mean(data_LEGs[['TAG_FourthG']])
LEGs_tag_c <- mean(data_LEGs[['TAG_FourthC']])

#Create dataframes - ALL
total.taa <- rbind(taa_a, taa_t, taa_g, taa_c)
colnames(total.taa) <- c('Frequency')
rownames(total.taa) <- c('A', 'T', 'G', 'C')

total.tga <- rbind(tga_a, tga_t, tga_g, tga_c)
colnames(total.tga) <- c('Frequency')
rownames(total.tga) <- c('A', 'T', 'G', 'C')

total.tag <- rbind(tag_a, tag_t, tag_g, tag_c)
colnames(total.tag) <- c('Frequency')
rownames(total.tag) <- c('A', 'T', 'G', 'C')

all <- cbind(total.taa, total.tga, total.tag)
colnames(all) <- c('TAA', 'TGA', 'TAG')

#Create dataframes - HEGs
hegs.total.taa <- rbind(HEGs_taa_a, HEGs_taa_t, HEGs_taa_g, HEGs_taa_c)
colnames(hegs.total.taa) <- c('Frequency')
rownames(hegs.total.taa) <- c('A', 'T', 'G', 'C')

hegs.total.tga <- rbind(HEGs_tga_a, HEGs_tga_t, HEGs_tga_g, HEGs_tga_c)
colnames(hegs.total.tga) <- c('Frequency')
rownames(hegs.total.tga) <- c('A', 'T', 'G', 'C')

hegs.total.tag <- rbind(HEGs_tag_a, HEGs_tag_t, HEGs_tag_g, HEGs_tag_c)
colnames(hegs.total.tag) <- c('Frequency')
rownames(hegs.total.tag) <- c('A', 'T', 'G', 'C')

hegs <- cbind(hegs.total.taa, hegs.total.tga, hegs.total.tag)
colnames(hegs) <- c('TAA', 'TGA', 'TAG')

#Create dataframes - LEGs
legs.total.taa <- rbind(LEGs_taa_a, LEGs_taa_t, LEGs_taa_g, LEGs_taa_c)
colnames(legs.total.taa) <- c('Frequency')
rownames(legs.total.taa) <- c('A', 'T', 'G', 'C')

legs.total.tga <- rbind(LEGs_tga_a, LEGs_tga_t, LEGs_tga_g, LEGs_tga_c)
colnames(legs.total.tga) <- c('Frequency')
rownames(legs.total.tga) <- c('A', 'T', 'G', 'C')

legs.total.tag <- rbind(LEGs_tag_a, LEGs_tag_t, LEGs_tag_g, LEGs_tag_c)
colnames(legs.total.tag) <- c('Frequency')
rownames(legs.total.tag) <- c('A', 'T', 'G', 'C')

legs <- cbind(legs.total.taa, legs.total.tga, legs.total.tag)
colnames(legs) <- c('TAA', 'TGA', 'TAG')

#Plot
par(mfrow=c(1,3), cex=1.5)

barplot(as.matrix((all)), beside=TRUE, col=c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), 
        border=NA, ylim=c(0, 0.6), ylab='Frequency of 4th Base', xlab='Primary Stop',
        main= 'All Genes', cex.lab =1.5, cex.axis=1.5, cex.names = 1.5, cex.main=2)
legend("topright", 
       legend = c("A", "T", "G", "C"),
       fill = c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), bty = "n", cex=1.2)

barplot(as.matrix((hegs)), beside=TRUE, col=c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), 
        border=NA, ylim=c(0, 0.6), ylab='Frequency of 4th Base', xlab='Primary Stop',
        main= 'HEGs', cex.lab =1.5, cex.axis=1.5, cex.names = 1.5, cex.main=2)
legend("topright", 
       legend = c("A", "T", "G", "C"),
       fill = c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), bty = "n", cex=1.2)

barplot(as.matrix((legs)), beside=TRUE, col=c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), 
        border=NA, ylim=c(0, 0.6), ylab='Frequency of 4th Base', xlab='Primary Stop',
        main= 'LEGs', cex.lab =1.5, cex.axis=1.5, cex.names = 1.5, cex.main=2)
legend("topright", 
       legend = c("A", "T", "G", "C"),
       fill = c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), bty = "n", cex=1.2)

#kruskall tests
k.all.taa <- kruskal.test(list(data_all[['TAA_FourthA']], data_all[['TAA_FourthT']], data_all[['TAA_FourthG']], data_all[['TAA_FourthC']]))
k.all.tga <- kruskal.test(list(data_all[['TGA_FourthA']], data_all[['TGA_FourthT']], data_all[['TGA_FourthG']], data_all[['TGA_FourthC']]))
k.all.tag <- kruskal.test(list(data_all[['TAG_FourthA']], data_all[['TAG_FourthT']], data_all[['TAG_FourthG']], data_all[['TAG_FourthC']]))

k.HEGs.taa <- kruskal.test(list(data_HEGs[['TAA_FourthA']], data_HEGs[['TAA_FourthT']], data_HEGs[['TAA_FourthG']], data_HEGs[['TAA_FourthC']]))
k.HEGs.tga <- kruskal.test(list(data_HEGs[['TGA_FourthA']], data_HEGs[['TGA_FourthT']], data_HEGs[['TGA_FourthG']], data_HEGs[['TGA_FourthC']]))
k.HEGs.tag <- kruskal.test(list(data_HEGs[['TAG_FourthA']], data_HEGs[['TAG_FourthT']], data_HEGs[['TAG_FourthG']], data_HEGs[['TAG_FourthC']]))

k.LEGs.taa <- kruskal.test(list(data_LEGs[['TAA_FourthA']], data_LEGs[['TAA_FourthT']], data_LEGs[['TAA_FourthG']], data_LEGs[['TAA_FourthC']]))
k.LEGs.tga <- kruskal.test(list(data_LEGs[['TGA_FourthA']], data_LEGs[['TGA_FourthT']], data_LEGs[['TGA_FourthG']], data_LEGs[['TGA_FourthC']]))
k.LEGs.tag <- kruskal.test(list(data_LEGs[['TAG_FourthA']], data_LEGs[['TAG_FourthT']], data_LEGs[['TAG_FourthG']], data_LEGs[['TAG_FourthC']]))

#Wilcoxons
w.all.tga <- wilcox.test(data_all[['TGA_FourthT']], data_all[['TGA_FourthA']], paired = TRUE, alternative = "greater")

w.HEGs.taa <- wilcox.test(data_HEGs[['TAA_FourthT']], data_HEGs[['TAA_FourthA']], paired = TRUE, alternative = "greater")
w.HEGs.tga <- wilcox.test(data_HEGs[['TGA_FourthT']], data_HEGs[['TGA_FourthA']], paired = TRUE, alternative = "greater")
w.HEGs.tag <- wilcox.test(data_HEGs[['TAG_FourthT']], data_HEGs[['TAG_FourthA']], paired = TRUE, alternative = "greater")

w.LEGs.tga <- wilcox.test(data_LEGs[['TGA_FourthT']], data_LEGs[['TGA_FourthA']], paired = TRUE, alternative = "greater")

g <- wilcox.test(data_all[['TAA_FourthT']], data_all[['TAG_FourthT']], alternative = "greater")
h <- wilcox.test(data_all[['TGA_FourthT']], data_all[['TAG_FourthT']], alternative = "greater")
