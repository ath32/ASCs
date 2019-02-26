### SOURCE

#Change source if needed (all genes, HEGs, LEGs)
source_all <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/10.4_T_codons_all.csv"
data_all = read.csv(source_all, header = TRUE)

source_tga <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/10.4_T_codons_tga.csv"
data_tga = read.csv(source_tga, header = TRUE)

source_taa <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/10.4_T_codons_taa.csv"
data_taa = read.csv(source_taa, header = TRUE)

source_tag <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/10.4_T_codons_tag.csv"
data_tag = read.csv(source_tag, header = TRUE)

source_hegs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/10.4_T_codons_hegs.csv"
data_hegs = read.csv(source_hegs, header = TRUE)

source_legs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/10.4_T_codons_legs.csv"
data_legs = read.csv(source_legs, header = TRUE)


#Create dataframe for each codon - accession, number, codon
#All
all_taa <- data_all[,c("Accession","taa_ratio")]
all_taa$Codon <- 'TAA'
colnames(all_taa) <- c('Accession', 'Score', 'Codon')

all_tga <- data_all[,c("Accession","tga_ratio")]
all_tga$Codon <- 'TGA'
colnames(all_tga) <- c('Accession', 'Score', 'Codon')

all_tag <- data_all[,c("Accession","tag_ratio")]
all_tag$Codon <- 'TAG'
colnames(all_tag) <- c('Accession', 'Score', 'Codon')

all_ttt <- data_all[,c("Accession","ttt_ratio")]
all_ttt$Codon <- 'TTT'
colnames(all_ttt) <- c('Accession', 'Score', 'Codon')

all_ttc <- data_all[,c("Accession","ttc_ratio")]
all_ttc$Codon <- 'TTC'
colnames(all_ttc) <- c('Accession', 'Score', 'Codon')

all_tta <- data_all[,c("Accession","tta_ratio")]
all_tta$Codon <- 'TTA'
colnames(all_tta) <- c('Accession', 'Score', 'Codon')

all_ttg <- data_all[,c("Accession","ttg_ratio")]
all_ttg$Codon <- 'TTG'
colnames(all_ttg) <- c('Accession', 'Score', 'Codon')

all_tct <- data_all[,c("Accession","tct_ratio")]
all_tct$Codon <- 'TCT'
colnames(all_tct) <- c('Accession', 'Score', 'Codon')

all_tcc <- data_all[,c("Accession","tcc_ratio")]
all_tcc$Codon <- 'TCC'
colnames(all_tcc) <- c('Accession', 'Score', 'Codon')

all_tca <- data_all[,c("Accession","tca_ratio")]
all_tca$Codon <- 'TCA'
colnames(all_tca) <- c('Accession', 'Score', 'Codon')

all_tcg <- data_all[,c("Accession","tcg_ratio")]
all_tcg$Codon <- 'TCG'
colnames(all_tcg) <- c('Accession', 'Score', 'Codon')

all_tat <- data_all[,c("Accession","tat_ratio")]
all_tat$Codon <- 'TAT'
colnames(all_tat) <- c('Accession', 'Score', 'Codon')

all_tac <- data_all[,c("Accession","tac_ratio")]
all_tac$Codon <- 'TAC'
colnames(all_tac) <- c('Accession', 'Score', 'Codon')

all_tgt <- data_all[,c("Accession","tgt_ratio")]
all_tgt$Codon <- 'TGT'
colnames(all_tgt) <- c('Accession', 'Score', 'Codon')

all_tgc <- data_all[,c("Accession","tgc_ratio")]
all_tgc$Codon <- 'TGC'
colnames(all_tgc) <- c('Accession', 'Score', 'Codon')

all_tgg <- data_all[,c("Accession","tgg_ratio")]
all_tgg$Codon <- 'TGG'
colnames(all_tgg) <- c('Accession', 'Score', 'Codon')

combined_all <- rbind(all_taa, all_tga, all_tag, all_ttt, all_ttc, all_tta, all_ttg, 
                      all_tct, all_tcc, all_tca, all_tcg, all_tat, all_tac, all_tgt,
                      all_tgc, all_tgg)

#TGA-terminating genes
tga_taa <- data_tga[,c("Accession","taa_ratio")]
tga_taa$Codon <- 'TAA'
colnames(tga_taa) <- c('Accession', 'Score', 'Codon')

tga_tga <- data_tga[,c("Accession","tga_ratio")]
tga_tga$Codon <- 'TGA'
colnames(tga_tga) <- c('Accession', 'Score', 'Codon')

tga_tag <- data_tga[,c("Accession","tag_ratio")]
tga_tag$Codon <- 'TAG'
colnames(tga_tag) <- c('Accession', 'Score', 'Codon')

tga_ttt <- data_tga[,c("Accession","ttt_ratio")]
tga_ttt$Codon <- 'TTT'
colnames(tga_ttt) <- c('Accession', 'Score', 'Codon')

tga_ttc <- data_tga[,c("Accession","ttc_ratio")]
tga_ttc$Codon <- 'TTC'
colnames(tga_ttc) <- c('Accession', 'Score', 'Codon')

tga_tta <- data_tga[,c("Accession","tta_ratio")]
tga_tta$Codon <- 'TTA'
colnames(tga_tta) <- c('Accession', 'Score', 'Codon')

tga_ttg <- data_tga[,c("Accession","ttg_ratio")]
tga_ttg$Codon <- 'TTG'
colnames(tga_ttg) <- c('Accession', 'Score', 'Codon')

tga_tct <- data_tga[,c("Accession","tct_ratio")]
tga_tct$Codon <- 'TCT'
colnames(tga_tct) <- c('Accession', 'Score', 'Codon')

tga_tcc <- data_tga[,c("Accession","tcc_ratio")]
tga_tcc$Codon <- 'TCC'
colnames(tga_tcc) <- c('Accession', 'Score', 'Codon')

tga_tca <- data_tga[,c("Accession","tca_ratio")]
tga_tca$Codon <- 'TCA'
colnames(tga_tca) <- c('Accession', 'Score', 'Codon')

tga_tcg <- data_tga[,c("Accession","tcg_ratio")]
tga_tcg$Codon <- 'TCG'
colnames(tga_tcg) <- c('Accession', 'Score', 'Codon')

tga_tat <- data_tga[,c("Accession","tat_ratio")]
tga_tat$Codon <- 'TAT'
colnames(tga_tat) <- c('Accession', 'Score', 'Codon')

tga_tac <- data_tga[,c("Accession","tac_ratio")]
tga_tac$Codon <- 'TAC'
colnames(tga_tac) <- c('Accession', 'Score', 'Codon')

tga_tgt <- data_tga[,c("Accession","tgt_ratio")]
tga_tgt$Codon <- 'TGT'
colnames(tga_tgt) <- c('Accession', 'Score', 'Codon')

tga_tgc <- data_tga[,c("Accession","tgc_ratio")]
tga_tgc$Codon <- 'TGC'
colnames(tga_tgc) <- c('Accession', 'Score', 'Codon')

tga_tgg <- data_tga[,c("Accession","tgg_ratio")]
tga_tgg$Codon <- 'TGG'
colnames(tga_tgg) <- c('Accession', 'Score', 'Codon')

combined_tga <- rbind(tga_taa, tga_tga, tga_tag, tga_ttt, tga_ttc, tga_tta, tga_ttg, 
                      tga_tct, tga_tcc, tga_tca, tga_tcg, tga_tat, tga_tac, tga_tgt,
                      tga_tgc, tga_tgg)

#TAA-terminating genes
taa_taa <- data_taa[,c("Accession","taa_ratio")]
taa_taa$Codon <- 'TAA'
colnames(taa_taa) <- c('Accession', 'Score', 'Codon')

taa_tga <- data_taa[,c("Accession","tga_ratio")]
taa_tga$Codon <- 'TGA'
colnames(taa_tga) <- c('Accession', 'Score', 'Codon')

taa_tag <- data_taa[,c("Accession","tag_ratio")]
taa_tag$Codon <- 'TAG'
colnames(taa_tag) <- c('Accession', 'Score', 'Codon')

taa_ttt <- data_taa[,c("Accession","ttt_ratio")]
taa_ttt$Codon <- 'TTT'
colnames(taa_ttt) <- c('Accession', 'Score', 'Codon')

taa_ttc <- data_taa[,c("Accession","ttc_ratio")]
taa_ttc$Codon <- 'TTC'
colnames(taa_ttc) <- c('Accession', 'Score', 'Codon')

taa_tta <- data_taa[,c("Accession","tta_ratio")]
taa_tta$Codon <- 'TTA'
colnames(taa_tta) <- c('Accession', 'Score', 'Codon')

taa_ttg <- data_taa[,c("Accession","ttg_ratio")]
taa_ttg$Codon <- 'TTG'
colnames(taa_ttg) <- c('Accession', 'Score', 'Codon')

taa_tct <- data_taa[,c("Accession","tct_ratio")]
taa_tct$Codon <- 'TCT'
colnames(taa_tct) <- c('Accession', 'Score', 'Codon')

taa_tcc <- data_taa[,c("Accession","tcc_ratio")]
taa_tcc$Codon <- 'TCC'
colnames(taa_tcc) <- c('Accession', 'Score', 'Codon')

taa_tca <- data_taa[,c("Accession","tca_ratio")]
taa_tca$Codon <- 'TCA'
colnames(taa_tca) <- c('Accession', 'Score', 'Codon')

taa_tcg <- data_taa[,c("Accession","tcg_ratio")]
taa_tcg$Codon <- 'TCG'
colnames(taa_tcg) <- c('Accession', 'Score', 'Codon')

taa_tat <- data_taa[,c("Accession","tat_ratio")]
taa_tat$Codon <- 'TAT'
colnames(taa_tat) <- c('Accession', 'Score', 'Codon')

taa_tac <- data_taa[,c("Accession","tac_ratio")]
taa_tac$Codon <- 'TAC'
colnames(taa_tac) <- c('Accession', 'Score', 'Codon')

taa_tgt <- data_taa[,c("Accession","tgt_ratio")]
taa_tgt$Codon <- 'TGT'
colnames(taa_tgt) <- c('Accession', 'Score', 'Codon')

taa_tgc <- data_taa[,c("Accession","tgc_ratio")]
taa_tgc$Codon <- 'TGC'
colnames(taa_tgc) <- c('Accession', 'Score', 'Codon')

taa_tgg <- data_taa[,c("Accession","tgg_ratio")]
taa_tgg$Codon <- 'TGG'
colnames(taa_tgg) <- c('Accession', 'Score', 'Codon')

combined_taa <- rbind(taa_taa, taa_tga, taa_tag, taa_ttt, taa_ttc, taa_tta, taa_ttg, 
                      taa_tct, taa_tcc, taa_tca, taa_tcg, taa_tat, taa_tac, taa_tgt,
                      taa_tgc, taa_tgg)

#TAG-terminating genes
tag_taa <- data_tag[,c("Accession","taa_ratio")]
tag_taa$Codon <- 'TAA'
colnames(tag_taa) <- c('Accession', 'Score', 'Codon')

tag_tga <- data_tag[,c("Accession","tga_ratio")]
tag_tga$Codon <- 'TGA'
colnames(tag_tga) <- c('Accession', 'Score', 'Codon')

tag_tag <- data_tag[,c("Accession","tag_ratio")]
tag_tag$Codon <- 'TAG'
colnames(tag_tag) <- c('Accession', 'Score', 'Codon')

tag_ttt <- data_tag[,c("Accession","ttt_ratio")]
tag_ttt$Codon <- 'TTT'
colnames(tag_ttt) <- c('Accession', 'Score', 'Codon')

tag_ttc <- data_tag[,c("Accession","ttc_ratio")]
tag_ttc$Codon <- 'TTC'
colnames(tag_ttc) <- c('Accession', 'Score', 'Codon')

tag_tta <- data_tag[,c("Accession","tta_ratio")]
tag_tta$Codon <- 'TTA'
colnames(tag_tta) <- c('Accession', 'Score', 'Codon')

tag_ttg <- data_tag[,c("Accession","ttg_ratio")]
tag_ttg$Codon <- 'TTG'
colnames(tag_ttg) <- c('Accession', 'Score', 'Codon')

tag_tct <- data_tag[,c("Accession","tct_ratio")]
tag_tct$Codon <- 'TCT'
colnames(tag_tct) <- c('Accession', 'Score', 'Codon')

tag_tcc <- data_tag[,c("Accession","tcc_ratio")]
tag_tcc$Codon <- 'TCC'
colnames(tag_tcc) <- c('Accession', 'Score', 'Codon')

tag_tca <- data_tag[,c("Accession","tca_ratio")]
tag_tca$Codon <- 'TCA'
colnames(tag_tca) <- c('Accession', 'Score', 'Codon')

tag_tcg <- data_tag[,c("Accession","tcg_ratio")]
tag_tcg$Codon <- 'TCG'
colnames(tag_tcg) <- c('Accession', 'Score', 'Codon')

tag_tat <- data_tag[,c("Accession","tat_ratio")]
tag_tat$Codon <- 'TAT'
colnames(tag_tat) <- c('Accession', 'Score', 'Codon')

tag_tac <- data_tag[,c("Accession","tac_ratio")]
tag_tac$Codon <- 'TAC'
colnames(tag_tac) <- c('Accession', 'Score', 'Codon')

tag_tgt <- data_tag[,c("Accession","tgt_ratio")]
tag_tgt$Codon <- 'TGT'
colnames(tag_tgt) <- c('Accession', 'Score', 'Codon')

tag_tgc <- data_tag[,c("Accession","tgc_ratio")]
tag_tgc$Codon <- 'TGC'
colnames(tag_tgc) <- c('Accession', 'Score', 'Codon')

tag_tgg <- data_tag[,c("Accession","tgg_ratio")]
tag_tgg$Codon <- 'TGG'
colnames(tag_tgg) <- c('Accession', 'Score', 'Codon')

combined_tag <- rbind(tag_taa, tag_tga, tag_tag, tag_ttt, tag_ttc, tag_tta, tag_ttg, 
                      tag_tct, tag_tcc, tag_tca, tag_tcg, tag_tat, tag_tac, tag_tgt,
                      tag_tgc, tag_tgg)

#LEGs
legs_taa <- data_legs[,c("Accession","taa_ratio")]
legs_taa$Codon <- 'TAA'
colnames(legs_taa) <- c('Accession', 'Score', 'Codon')

legs_tga <- data_legs[,c("Accession","tga_ratio")]
legs_tga$Codon <- 'TGA'
colnames(legs_tga) <- c('Accession', 'Score', 'Codon')

legs_tag <- data_legs[,c("Accession","tag_ratio")]
legs_tag$Codon <- 'TAG'
colnames(legs_tag) <- c('Accession', 'Score', 'Codon')

legs_ttt <- data_legs[,c("Accession","ttt_ratio")]
legs_ttt$Codon <- 'TTT'
colnames(legs_ttt) <- c('Accession', 'Score', 'Codon')

legs_ttc <- data_legs[,c("Accession","ttc_ratio")]
legs_ttc$Codon <- 'TTC'
colnames(legs_ttc) <- c('Accession', 'Score', 'Codon')

legs_tta <- data_legs[,c("Accession","tta_ratio")]
legs_tta$Codon <- 'TTA'
colnames(legs_tta) <- c('Accession', 'Score', 'Codon')

legs_ttg <- data_legs[,c("Accession","ttg_ratio")]
legs_ttg$Codon <- 'TTG'
colnames(legs_ttg) <- c('Accession', 'Score', 'Codon')

legs_tct <- data_legs[,c("Accession","tct_ratio")]
legs_tct$Codon <- 'TCT'
colnames(legs_tct) <- c('Accession', 'Score', 'Codon')

legs_tcc <- data_legs[,c("Accession","tcc_ratio")]
legs_tcc$Codon <- 'TCC'
colnames(legs_tcc) <- c('Accession', 'Score', 'Codon')

legs_tca <- data_legs[,c("Accession","tca_ratio")]
legs_tca$Codon <- 'TCA'
colnames(legs_tca) <- c('Accession', 'Score', 'Codon')

legs_tcg <- data_legs[,c("Accession","tcg_ratio")]
legs_tcg$Codon <- 'TCG'
colnames(legs_tcg) <- c('Accession', 'Score', 'Codon')

legs_tat <- data_legs[,c("Accession","tat_ratio")]
legs_tat$Codon <- 'TAT'
colnames(legs_tat) <- c('Accession', 'Score', 'Codon')

legs_tac <- data_legs[,c("Accession","tac_ratio")]
legs_tac$Codon <- 'TAC'
colnames(legs_tac) <- c('Accession', 'Score', 'Codon')

legs_tgt <- data_legs[,c("Accession","tgt_ratio")]
legs_tgt$Codon <- 'TGT'
colnames(legs_tgt) <- c('Accession', 'Score', 'Codon')

legs_tgc <- data_legs[,c("Accession","tgc_ratio")]
legs_tgc$Codon <- 'TGC'
colnames(legs_tgc) <- c('Accession', 'Score', 'Codon')

legs_tgg <- data_legs[,c("Accession","tgg_ratio")]
legs_tgg$Codon <- 'TGG'
colnames(legs_tgg) <- c('Accession', 'Score', 'Codon')

combined_legs <- rbind(legs_taa, legs_tga, legs_tag, legs_ttt, legs_ttc, legs_tta, legs_ttg, 
                       legs_tct, legs_tcc, legs_tca, legs_tcg, legs_tat, legs_tac, legs_tgt,
                       legs_tgc, legs_tgg)

#HEGS
hegs_taa <- data_hegs[,c("Accession","taa_ratio")]
hegs_taa$Codon <- 'TAA'
colnames(hegs_taa) <- c('Accession', 'Score', 'Codon')

hegs_tga <- data_hegs[,c("Accession","tga_ratio")]
hegs_tga$Codon <- 'TGA'
colnames(hegs_tga) <- c('Accession', 'Score', 'Codon')

hegs_tag <- data_hegs[,c("Accession","tag_ratio")]
hegs_tag$Codon <- 'TAG'
colnames(hegs_tag) <- c('Accession', 'Score', 'Codon')

hegs_ttt <- data_hegs[,c("Accession","ttt_ratio")]
hegs_ttt$Codon <- 'TTT'
colnames(hegs_ttt) <- c('Accession', 'Score', 'Codon')

hegs_ttc <- data_hegs[,c("Accession","ttc_ratio")]
hegs_ttc$Codon <- 'TTC'
colnames(hegs_ttc) <- c('Accession', 'Score', 'Codon')

hegs_tta <- data_hegs[,c("Accession","tta_ratio")]
hegs_tta$Codon <- 'TTA'
colnames(hegs_tta) <- c('Accession', 'Score', 'Codon')

hegs_ttg <- data_hegs[,c("Accession","ttg_ratio")]
hegs_ttg$Codon <- 'TTG'
colnames(hegs_ttg) <- c('Accession', 'Score', 'Codon')

hegs_tct <- data_hegs[,c("Accession","tct_ratio")]
hegs_tct$Codon <- 'TCT'
colnames(hegs_tct) <- c('Accession', 'Score', 'Codon')

hegs_tcc <- data_hegs[,c("Accession","tcc_ratio")]
hegs_tcc$Codon <- 'TCC'
colnames(hegs_tcc) <- c('Accession', 'Score', 'Codon')

hegs_tca <- data_hegs[,c("Accession","tca_ratio")]
hegs_tca$Codon <- 'TCA'
colnames(hegs_tca) <- c('Accession', 'Score', 'Codon')

hegs_tcg <- data_hegs[,c("Accession","tcg_ratio")]
hegs_tcg$Codon <- 'TCG'
colnames(hegs_tcg) <- c('Accession', 'Score', 'Codon')

hegs_tat <- data_hegs[,c("Accession","tat_ratio")]
hegs_tat$Codon <- 'TAT'
colnames(hegs_tat) <- c('Accession', 'Score', 'Codon')

hegs_tac <- data_hegs[,c("Accession","tac_ratio")]
hegs_tac$Codon <- 'TAC'
colnames(hegs_tac) <- c('Accession', 'Score', 'Codon')

hegs_tgt <- data_hegs[,c("Accession","tgt_ratio")]
hegs_tgt$Codon <- 'TGT'
colnames(hegs_tgt) <- c('Accession', 'Score', 'Codon')

hegs_tgc <- data_hegs[,c("Accession","tgc_ratio")]
hegs_tgc$Codon <- 'TGC'
colnames(hegs_tgc) <- c('Accession', 'Score', 'Codon')

hegs_tgg <- data_hegs[,c("Accession","tgg_ratio")]
hegs_tgg$Codon <- 'TGG'
colnames(hegs_tgg) <- c('Accession', 'Score', 'Codon')

combined_hegs <- rbind(hegs_taa, hegs_tga, hegs_tag, hegs_ttt, hegs_ttc, hegs_tta, hegs_ttg, 
                       hegs_tct, hegs_tcc, hegs_tca, hegs_tcg, hegs_tat, hegs_tac, hegs_tgt,
                       hegs_tgc, hegs_tgg)

#PLOTS

par(mfrow=c(2,3), cex=1.1, mai = c(0.6, 1, 1, 1))

#all
boxplot(combined_all$Score ~ combined_all$Codon, main='All genes', outline=FALSE, las=2, ylab='Enrichment score', ylim=c(-1,3), xlab = 'Codon', col=c('lightblue', 'white', 'lightblue', 'white', 'white', 
                                                                                                                                             'white', 'white', 'white', 'lightblue', 'white', 'white',
                                                                                                                                             'white', 'white', 'white', 'white', 'white'))
abline(h=0, lty=2, lwd=2, col='red')

#hegs
boxplot(combined_hegs$Score ~ combined_hegs$Codon, main='HEGs',  outline=FALSE, las=2, ylab='Enrichment score', xlab = 'Codon', col=c('lightblue', 'white', 'lightblue', 'white', 'white', 
                                                                                                                                      'white', 'white', 'white', 'lightblue', 'white', 'white',
                                                                                                                                      'white', 'white', 'white', 'white', 'white'))
abline(h=0, lty=2, lwd=2, col='red')

#legs
boxplot(combined_legs$Score ~ combined_legs$Codon, main='LEGs',  outline=FALSE, las=2, ylab='Enrichment score', xlab = 'Codon', col=c('lightblue', 'white', 'lightblue', 'white', 'white', 
                                                                                                                                'white', 'white', 'white', 'lightblue', 'white', 'white',
                                                                                                                                'white', 'white', 'white', 'white', 'white'))
abline(h=0, lty=2, lwd=2, col='red')

#TAA
boxplot(combined_taa$Score ~ combined_taa$Codon, main='TAA-terminating genes', las=2, outline=FALSE, ylab='Enrichment score', xlab = 'Codon', col=c('lightblue', 'white', 'lightblue', 'white', 'white', 
                                                                                                                                             'white', 'white', 'white', 'lightblue', 'white', 'white',
                                                                                                                                             'white', 'white', 'white', 'white', 'white'))
abline(h=0, lty=2, lwd=2, col='red')

#TGA
boxplot(combined_tga$Score ~ combined_tga$Codon, main='TGA-terminating genes', las=2, outline=FALSE, ylab='Enrichment score', xlab = 'Codon', col=c('lightblue', 'white', 'lightblue', 'white', 'white', 
                                                                                                                                                    'white', 'white', 'white', 'lightblue', 'white', 'white',
                                                                                                                                                    'white', 'white', 'white', 'white', 'white'))
abline(h=0, lty=2, lwd=2, col='red')

#TAG
boxplot(combined_tag$Score ~ combined_tag$Codon, main='TAG-terminating genes', las=2, outline=FALSE, ylab='Enrichment score', xlab = 'Codon', col=c('lightblue', 'white', 'lightblue', 'white', 'white', 
                                                                                                                                                    'white', 'white', 'white', 'lightblue', 'white', 'white',
                                                                                                                                                    'white', 'white', 'white', 'white', 'white'))
abline(h=0, lty=2, lwd=2, col='red')


