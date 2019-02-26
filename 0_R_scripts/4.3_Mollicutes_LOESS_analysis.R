### SOURCE ### 

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.5_Frequencies_all_2.csv"
data = read.csv(source, header = TRUE)

observed_4 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.5_Frequencies_mollicutes_2.csv"
data_4 = read.csv(observed_4, header = TRUE)

#Plot LOESS models at each positions

par(mfrow=c(2,4), cex=1.3, cex.main=2)

scatter.smooth(data[['GC3']], data[['tga_position_0']], pch=16, col="#1b9e77", xlab = "GC3 content (%)", ylab = "Codon frequency", main='Position +0')
loess.smooth(data[['GC3']], data[['tga_position_0']], span = 2/3, degree = 1, col="#1b9e77")

scatter.smooth(data[['GC3']], data[['tga_position_1']], pch=16, col="#1b9e77", xlab = "GC3 content (%)", ylab = "Codon frequency", main='Position +1')
loess.smooth(data[['GC3']], data[['tga_position_1']], span = 2/3, degree = 1, col="#1b9e77")

scatter.smooth(data[['GC3']], data[['tga_position_2']], pch=16, col="#d95f02", xlab = "GC3 content (%)", ylab = "Codon frequency", main='Position +2')
loess.smooth(data[['GC3']], data[['tga_position_2']], span = 2/3, degree = 1, col="#d95f02")

scatter.smooth(data[['GC3']], data[['tga_position_3']], pch=16, col="#7570b3", xlab = "GC3 content (%)", ylab = "Codon frequency", main='Position +3')
loess.smooth(data[['GC3']], data[['tga_position_3']], span = 2/3, degree = 1, col="#d7570b3")

scatter.smooth(data[['GC3']], data[['tga_position_4']], pch=16, col="#e7298a", xlab = "GC3 content (%)", ylab = "Codon frequency", main='Position +4')
loess.smooth(data[['GC3']], data[['tga_position_4']], span = 2/3, degree = 1, col="#e7298a")

scatter.smooth(data[['GC3']], data[['tga_position_5']], pch=16, col="#66a61e", xlab = "GC3 content (%)", ylab = "Codon frequency", main='Position +5')
loess.smooth(data[['GC3']], data[['tga_position_5']], span = 2/3, degree = 1, col="#66a61e")

scatter.smooth(data[['GC3']], data[['tga_position_6']], pch=16, col="#e6ab02", xlab = "GC3 content (%)", ylab = "Codon frequency", main='Position +6')
loess.smooth(data[['GC3']], data[['tga_position_6']], span = 2/3, degree = 1, col="#e6ab02")

#Predictions

#Get IDs
colnames(data)
ids = colnames(data)[c(2:length(colnames(data))-2)]
newids = ids[c(-1)]
gc3 = data$GC3

ids_4 = colnames(data_4)[2:length(colnames(data))-2]
newids_4 = ids_4[c(-1)]
gc3_4 = data_4$GC3

#Get models
models <- list()

for (id in newids) {
  i <- data[[id]]
  gc3 <- data[['GC3']]
  models[[id]] <- loess(i ~ gc3, span=2/3, control=loess.control(surface="direct"))
} 

#Make predictions
gc3_vals = c()

for (i in 1:nrow(data_4)) {
  gc3_val <- data_4[i,"GC3"]
  gc3_vals <- c(gc3_vals, gc3_val)
}   

predictions <- list()

for (id in newids_4) {
  predictions[[id]] <- predict(object = models[[id]], newdata = gc3_vals)
}

#Wilcoxons

results = list()

for (q in newids) {

  observed <- data_4[[q]]
  predicted <- predictions[[q]]

  wilcoxon <- wilcox.test(observed, predicted, paired = TRUE, alternative = "l")
  results[[q]] <- wilcoxon$p.value

}

results

#Output results as dataframe
df <- data.frame(matrix(unlist(results), nrow=64, byrow=T))
rownames(df) <- (c('ttt', 'ttc', 'tta', 'ttg', 'tct', 'tcc', 'tca', 'tcg', 'tat', 'tac', 'taa', 'tag', 'tgt', 'tgc', 'tga', 'tgg',
           'ctt', 'ctc', 'cta', 'ctg', 'cct', 'ccc', 'cca', 'ccg', 'cat', 'cac', 'caa', 'cag', 'cgt', 'cgc', 'cga', 'cgg',
           'att', 'atc', 'ata', 'atg', 'act', 'acc', 'aca', 'acg', 'aat', 'aac', 'aaa', 'aag', 'agt', 'agc', 'aga', 'agg',
           'gtt', 'gtc', 'gta', 'gtg', 'gct', 'gcc', 'gca', 'gcg', 'gat', 'gac', 'gaa', 'gag', 'ggt', 'ggc', 'gga', 'ggg'))
colnames(df) <- c('+0', '+1', '+2', '+3', '+4', '+5', '+6')

#Ordered dataframes (by position)
zero <- df[order(df[['+0']]),] 
one <- df[order(df[['+1']]),] 
two <- df[order(df[['+2']]),] 
three <- df[order(df[['+3']]),] 
four <- df[order(df[['+4']]),] 
five <- df[order(df[['+5']]),] 
six <- df[order(df[['+6']]),] 

