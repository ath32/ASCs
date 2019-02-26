### DATA PRESENTATION 
### To print, unhash the plot of interest

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/4.3_GC_matched.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/4.3_GC_matched_median.csv"
data = read.csv(source, header = TRUE)

#Create columns for differences
data['Diff1'] <- data['MF_1'] - data['PF_1']
data['Diff2'] <- data['MF_2'] - data['PF_2']
data['Diff3'] <- data['MF_3'] - data['PF_3']
data['Diff4'] <- data['MF_4'] - data['PF_4']
data['Diff5'] <- data['MF_5'] - data['PF_5']
data['Diff6'] <- data['MF_6'] - data['PF_6']

#Plots

par(mfrow=c(2,3), cex=1.3)

barplot(data[['Diff1']], ylab = 'GC matched (Mollicute - TT11)', ylim=c(-0.02, 0.06), main='Pos +1', col='lightblue')
barplot(data[['Diff2']], ylab = 'GC matched (Mollicute - TT11)', ylim=c(-0.02, 0.04), main='Pos +2', col='lightblue')
barplot(data[['Diff3']], ylab = 'GC matched (Mollicute - TT11)', ylim=c(-0.02, 0.04), main='Pos +3', col='lightblue')
barplot(data[['Diff4']], ylab = 'GC matched (Mollicute - TT11)', ylim=c(-0.02, 0.04), main='Pos +4', col='lightblue')
barplot(data[['Diff5']], ylab = 'GC matched (Mollicute - TT11)', ylim=c(-0.02, 0.1), main='Pos +5', col='lightblue')
barplot(data[['Diff6']], ylab = 'GC matched (Mollicute - TT11)', ylim=c(-0.02, 0.04), main='Pos +6', col='lightblue')

#Scatter
# plot(data[['GC3']], data[['Diff1']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
# abline(h=0, lty=2, col='red', lwd=3)
# title("Position +1")
# 
# plot(data[['GC3']], data[['Diff2']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
# abline(h=0, lty=2, col='red', lwd=3)
# title("Position +2")
# 
# plot(data[['GC3']], data[['Diff3']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
# abline(h=0, lty=2, col='red', lwd=3)
# title("Position +3")
# 
# plot(data[['GC3']], data[['Diff4']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
# abline(h=0, lty=2, col='red', lwd=3)
# title("Position +4")
# 
# plot(data[['GC3']], data[['Diff5']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
# abline(h=0, lty=2, col='red', lwd=3)
# title("Position +5")
# 
# plot(data[['GC3']], data[['Diff6']], pch=16, col="lightblue", xlab = "GC3 content (%)", ylab = "GC-matched difference (Mollicute - TT11)")
# abline(h=0, lty=2, col='red', lwd=3)
# title("Position +6")

#Statistics #TT4
w1 = wilcox.test(data[['MF_1']], data[['PF_1']], paired = TRUE, alternative = "l")
w2 = wilcox.test(data[['MF_2']], data[['PF_2']], paired = TRUE, alternative = "l")
w3 = wilcox.test(data[['MF_3']], data[['PF_3']], paired = TRUE, alternative = "l")
w4 = wilcox.test(data[['MF_4']], data[['PF_4']], paired = TRUE, alternative = "l")
w5 = wilcox.test(data[['MF_5']], data[['PF_5']], paired = TRUE, alternative = "l")
w6 = wilcox.test(data[['MF_6']], data[['PF_6']], paired = TRUE, alternative = "l")

w1
w2
w3
w4
w5
w6
