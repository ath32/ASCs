### SOURCE ###

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.5_Frequencies_all_positions.csv"
data = read.csv(source, header = TRUE)

observed_4 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/3.5_Frequencies_all_positions_mollicutes.csv"
data_4 = read.csv(observed_4, header = TRUE)

#Plotting LOESS for visualisation purposes - change codon to whatever you like

par(mfrow=c(1,1), cex=1.3, cex.main=2)

scatter.smooth(data[['GC3']], data[['ttt']], pch=16, col="#1b9e77", xlab = "GC3 content (%)", ylab = "Codon frequency", main='Position +0')
loess.smooth(data[['GC3']], data[['ttt']], span = 2/3, degree = 1, col="#1b9e77")

### IGNORE BELOW... SEE 4.3 FOR LOESS ANALYSIS ###

#Get IDs
colnames(data)
ids = colnames(data)[c(2:length(colnames(data))-2)]
newids = ids[c(-1)]
gc3 = data$GC3

#Get models
models <- list()

for (id in newids) {
  i <- data[[id]]
  gc3 <- data[['GC3']]
  models[[id]] <- loess(i ~ gc3, span=2/3, control=loess.control(surface="direct"))
} 


#Make predictions
gc3_vals = c()

for (j in 1:nrow(data_4)) {
  gc3_val <- data_4[j,"GC3"]
  gc3_vals <- c(gc3_vals, gc3_val)
}   

predictions <- list()

for (id in newids) {
  predictions[[id]] <- predict(object = models[[id]], newdata = gc3_vals)
}

#Wilcoxons

results = list()

for (id in newids) {
  
  observed <- data_4[[id]]
  predicted <- predictions[[id]]
  
  wilcoxon <- wilcox.test(observed, predicted, paired = TRUE, alternative = "l")
  results[[id]] <- wilcoxon$p.value
  
}

#Create dataframe
df <- data.frame(matrix(unlist(results), nrow=64, byrow=T))
rownames(df) <- (c('ttt', 'ttc', 'tta', 'ttg', 'tct', 'tcc', 'tca', 'tcg', 'tat', 'tac', 'taa', 'tag', 'tgt', 'tgc', 'tga', 'tgg',
                   'ctt', 'ctc', 'cta', 'ctg', 'cct', 'ccc', 'cca', 'ccg', 'cat', 'cac', 'caa', 'cag', 'cgt', 'cgc', 'cga', 'cgg',
                   'att', 'atc', 'ata', 'atg', 'act', 'acc', 'aca', 'acg', 'aat', 'aac', 'aaa', 'aag', 'agt', 'agc', 'aga', 'agg',
                   'gtt', 'gtc', 'gta', 'gtg', 'gct', 'gcc', 'gca', 'gcg', 'gat', 'gac', 'gaa', 'gag', 'ggt', 'ggc', 'gga', 'ggg'))
colnames(df) <- c('Frequency')

#Ordered dataframes
ordered <- df[order(df[['Frequency']]),] 
