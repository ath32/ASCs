### SOURCE ###

source_all <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.5_fifth_site_all.csv"
data_all = read.csv(source_all, header = TRUE)

source_HEGs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.5_fifth_site_hegs.csv"
data_HEGs = read.csv(source_HEGs, header = TRUE)

source_LEGs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.5_fifth_site_legs.csv"
data_LEGs = read.csv(source_LEGs, header = TRUE)

#Calculate means - ALL
taa_a <- mean(data_all[['TAA_fifthA']])
taa_t <- mean(data_all[['TAA_fifthT']])
taa_g <- mean(data_all[['TAA_fifthG']])
taa_c <- mean(data_all[['TAA_fifthC']])
tga_a <- mean(data_all[['TGA_fifthA']])
tga_t <- mean(data_all[['TGA_fifthT']])
tga_g <- mean(data_all[['TGA_fifthG']])
tga_c <- mean(data_all[['TGA_fifthC']])
tag_a <- mean(data_all[['TAG_fifthA']])
tag_t <- mean(data_all[['TAG_fifthT']])
tag_g <- mean(data_all[['TAG_fifthG']])
tag_c <- mean(data_all[['TAG_fifthC']])

#Calculate means - HEGs
HEGs_taa_a <- mean(data_HEGs[['TAA_fifthA']])
HEGs_taa_t <- mean(data_HEGs[['TAA_fifthT']])
HEGs_taa_g <- mean(data_HEGs[['TAA_fifthG']])
HEGs_taa_c <- mean(data_HEGs[['TAA_fifthC']])
HEGs_tga_a <- mean(data_HEGs[['TGA_fifthA']])
HEGs_tga_t <- mean(data_HEGs[['TGA_fifthT']])
HEGs_tga_g <- mean(data_HEGs[['TGA_fifthG']])
HEGs_tga_c <- mean(data_HEGs[['TGA_fifthC']])
HEGs_tag_a <- mean(data_HEGs[['TAG_fifthA']])
HEGs_tag_t <- mean(data_HEGs[['TAG_fifthT']])
HEGs_tag_g <- mean(data_HEGs[['TAG_fifthG']])
HEGs_tag_c <- mean(data_HEGs[['TAG_fifthC']])

#Calculate means - LEGs
LEGs_taa_a <- mean(data_LEGs[['TAA_fifthA']])
LEGs_taa_t <- mean(data_LEGs[['TAA_fifthT']])
LEGs_taa_g <- mean(data_LEGs[['TAA_fifthG']])
LEGs_taa_c <- mean(data_LEGs[['TAA_fifthC']])
LEGs_tga_a <- mean(data_LEGs[['TGA_fifthA']])
LEGs_tga_t <- mean(data_LEGs[['TGA_fifthT']])
LEGs_tga_g <- mean(data_LEGs[['TGA_fifthG']])
LEGs_tga_c <- mean(data_LEGs[['TGA_fifthC']])
LEGs_tag_a <- mean(data_LEGs[['TAG_fifthA']])
LEGs_tag_t <- mean(data_LEGs[['TAG_fifthT']])
LEGs_tag_g <- mean(data_LEGs[['TAG_fifthG']])
LEGs_tag_c <- mean(data_LEGs[['TAG_fifthC']])

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
par(mfrow=c(1,3), cex=1)

barplot(as.matrix((all)), beside=TRUE, col=c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), 
        border=NA, ylim=c(0, 0.6), ylab='Frequency of 5th Base in +4T genes', xlab='Primary Stop',
        main= 'All Genes', cex.lab =1.2)
legend("topright", 
       legend = c("A", "T", "G", "C"),
       fill = c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), bty = "n")

barplot(as.matrix((hegs)), beside=TRUE, col=c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), 
        border=NA, ylim=c(0, 0.6), ylab='Frequency of 5th Base in +4T genes', xlab='Primary Stop',
        main= 'HEGs', cex.lab =1.2)
legend("topright", 
       legend = c("A", "T", "G", "C"),
       fill = c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), bty = "n")

barplot(as.matrix((legs)), beside=TRUE, col=c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), 
        border=NA, ylim=c(0, 0.6), ylab='Frequency of 5th Base in +4T genes', xlab='Primary Stop',
        main= 'LEGs', cex.lab =1.2)
legend("topright", 
       legend = c("A", "T", "G", "C"),
       fill = c("#E1F5FE", "#01579B", '#81D4FA', '#0288D1'), bty = "n")

#Wilcoxons
tga.all.t <- wilcox.test(data_all[['TGA_fifthT']], data_all[['TGA_fifthA']], paired = TRUE, alternative = "greater")
tga.all.c <- wilcox.test(data_all[['TGA_fifthC']], data_all[['TGA_fifthA']], paired = TRUE, alternative = "greater")
taa.all.t <- wilcox.test(data_all[['TAA_fifthT']], data_all[['TAA_fifthA']], paired = TRUE, alternative = "greater")
taa.all.c <- wilcox.test(data_all[['TAA_fifthC']], data_all[['TAA_fifthA']], paired = TRUE, alternative = "greater")
tag.all.t <- wilcox.test(data_all[['TAG_fifthT']], data_all[['TAG_fifthA']], paired = TRUE, alternative = "greater")
tag.all.c <- wilcox.test(data_all[['TAG_fifthC']], data_all[['TAG_fifthA']], paired = TRUE, alternative = "greater")

taa.hegs.t <- wilcox.test(data_HEGs[['TAA_fifthT']], data_HEGs[['TAA_fifthC']], paired = TRUE, alternative = "greater")
tga.hegs.t <- wilcox.test(data_HEGs[['TGA_fifthT']], data_HEGs[['TGA_fifthC']], paired = TRUE, alternative = "greater")
tag.hegs.t <- wilcox.test(data_HEGs[['TAG_fifthT']], data_HEGs[['TAG_fifthC']], paired = TRUE, alternative = "greater")

taa.LEGs.tga <- wilcox.test(data_LEGs[['TAA_fifthT']], data_LEGs[['TAA_fifthA']], paired = TRUE, alternative = "greater")
tga.LEGs.tga <- wilcox.test(data_LEGs[['TGA_fifthT']], data_LEGs[['TGA_fifthA']], paired = TRUE, alternative = "greater")
tag.LEGs.tga <- wilcox.test(data_LEGs[['TAG_fifthT']], data_LEGs[['TAG_fifthG']], paired = TRUE, alternative = "greater")

k.legs.taa <- kruskal.test(list(data_LEGs[['TAA_fifthA']], data_LEGs[['TAA_fifthT']], data_LEGs[['TAA_fifthG']], data_LEGs[['TAA_fifthC']]))
k.legs.tga <- kruskal.test(list(data_LEGs[['TGA_fifthA']], data_LEGs[['TGA_fifthT']], data_LEGs[['TGA_fifthG']], data_LEGs[['TGA_fifthC']]))
k.legs.tag <- kruskal.test(list(data_LEGs[['TAG_fifthA']], data_LEGs[['TAG_fifthT']], data_LEGs[['TAG_fifthG']], data_LEGs[['TAG_fifthC']]))
