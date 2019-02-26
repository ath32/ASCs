### SOURCE ### Unhash source file of interest

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/2.2_Simulations.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_allgenes.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_HEGs.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/9.1_+4T_sims_LEGs.csv"
data = read.csv(source, header = TRUE)

#Zs
z.1 <- data[['P1_Zscore']]
z.2 <- data[['P2_Zscore']]
z.3 <- data[['P3_Zscore']]
z.4 <- data[['P4_Zscore']]
z.5 <- data[['P5_Zscore']]
z.6 <- data[['P6_Zscore']]

#Positive counts
p.1 <- length(z.1[which(z.1 > 0)])
p.2 <- length(z.2[which(z.2 > 0)])
p.3 <- length(z.3[which(z.3 > 0)])
p.4 <- length(z.4[which(z.4 > 0)])
p.5 <- length(z.5[which(z.5 > 0)])
p.6 <- length(z.6[which(z.6 > 0)])

#Negative counts
n.1 <- length(z.1[which(z.1 < 0)])
n.2 <- length(z.2[which(z.2 < 0)])
n.3 <- length(z.3[which(z.3 < 0)])
n.4 <- length(z.4[which(z.4 < 0)])
n.5 <- length(z.5[which(z.5 < 0)])
n.6 <- length(z.6[which(z.6 < 0)])

#Positive Binomials
pb.1 <- binom.test(p.1, length(z.1), p = 0.5, alternative = c("two.sided"))
pb.2 <- binom.test(p.2, length(z.2), p = 0.5, alternative = c("two.sided"))
pb.3 <- binom.test(p.3, length(z.3), p = 0.5, alternative = c("two.sided"))
pb.4 <- binom.test(p.4, length(z.4), p = 0.5, alternative = c("two.sided"))
pb.5 <- binom.test(p.5, length(z.5), p = 0.5, alternative = c("two.sided"))
pb.6 <- binom.test(p.6, length(z.6), p = 0.5, alternative = c("two.sided"))

#Negative Binomials
nb.1 <- binom.test(n.1, length(z.1), p = 0.5, alternative = c("two.sided"))
nb.2 <- binom.test(n.2, length(z.2), p = 0.5, alternative = c("two.sided"))
nb.3 <- binom.test(n.3, length(z.3), p = 0.5, alternative = c("two.sided"))
nb.4 <- binom.test(n.4, length(z.4), p = 0.5, alternative = c("two.sided"))
nb.5 <- binom.test(n.5, length(z.5), p = 0.5, alternative = c("two.sided"))
nb.6 <- binom.test(n.6, length(z.6), p = 0.5, alternative = c("two.sided"))

#One-tailed Positive Binomials
opb.1 <- binom.test(p.1, length(z.1), p = 0.5, alternative = c("l"))
opb.2 <- binom.test(p.2, length(z.2), p = 0.5, alternative = c("l"))
opb.3 <- binom.test(p.3, length(z.3), p = 0.5, alternative = c("l"))
opb.4 <- binom.test(p.4, length(z.4), p = 0.5, alternative = c("l"))
opb.5 <- binom.test(p.5, length(z.5), p = 0.5, alternative = c("l"))
opb.6 <- binom.test(p.6, length(z.6), p = 0.5, alternative = c("l"))

#Significant binomials
sd.1 <- length(z.1[which(abs(z.1) > 1.96)])
sd.2 <- length(z.2[which(abs(z.2) > 1.96)])
sd.3 <- length(z.3[which(abs(z.3) > 1.96)])
sd.4 <- length(z.4[which(abs(z.4) > 1.96)])
sd.5 <- length(z.5[which(abs(z.5) > 1.96)])
sd.6 <- length(z.6[which(abs(z.6) > 1.96)])

sob.1 <- binom.test(sd.1, length(z.1), p = 0.05, alternative = c("two.sided"))
sob.2 <- binom.test(sd.2, length(z.2), p = 0.05, alternative = c("two.sided"))
sob.3 <- binom.test(sd.3, length(z.3), p = 0.05, alternative = c("two.sided"))
sob.4 <- binom.test(sd.4, length(z.4), p = 0.05, alternative = c("two.sided"))
sob.5 <- binom.test(sd.5, length(z.5), p = 0.05, alternative = c("two.sided"))
sob.6 <- binom.test(sd.6, length(z.6), p = 0.05, alternative = c("two.sided"))

so.1 <- length(z.1[which(z.1 > 1.64)])
so.2 <- length(z.2[which(z.2 > 1.64)])
so.3 <- length(z.3[which(z.3 > 1.64)])
so.4 <- length(z.4[which(z.4 > 1.64)])
so.5 <- length(z.5[which(z.5 > 1.64)])
so.6 <- length(z.6[which(z.6 > 1.64)])

osob.1 <- binom.test(so.1, length(z.1), p = 0.05, alternative = c("l"))
osob.2 <- binom.test(so.2, length(z.2), p = 0.05, alternative = c("l"))
osob.3 <- binom.test(so.3, length(z.3), p = 0.05, alternative = c("l"))
osob.4 <- binom.test(so.4, length(z.4), p = 0.05, alternative = c("l"))
osob.5 <- binom.test(so.5, length(z.5), p = 0.05, alternative = c("l"))
osob.6 <- binom.test(so.6, length(z.6), p = 0.05, alternative = c("l"))

#Significant binomials - underenrichment
su.1 <- length(z.1[which(z.1 < -1.64)])
su.2 <- length(z.2[which(z.2 < -1.64)])
su.3 <- length(z.3[which(z.3 < -1.64)])
su.4 <- length(z.4[which(z.4 < -1.64)])
su.5 <- length(z.5[which(z.5 < -1.64)])
su.6 <- length(z.6[which(z.6 < -1.64)])

osub.1 <- binom.test(su.1, length(z.1), p = 0.05, alternative = c("g"))
osub.2 <- binom.test(su.2, length(z.2), p = 0.05, alternative = c("g"))
osub.3 <- binom.test(su.3, length(z.3), p = 0.05, alternative = c("g"))
osub.4 <- binom.test(su.4, length(z.4), p = 0.05, alternative = c("g"))
osub.5 <- binom.test(su.5, length(z.5), p = 0.05, alternative = c("g"))
osub.6 <- binom.test(su.6, length(z.6), p = 0.05, alternative = c("g"))
