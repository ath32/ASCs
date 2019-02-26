### SOURCE ###

source_total <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/6.3_Thirdstops_all.csv"
data_total = read.csv(source_total, header = TRUE)

source_HEGs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/6.3_Thirdstops_HEGs.csv"
data_HEGs = read.csv(source_HEGs, header = TRUE)

source_LEGs <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/6.3_Thirdstops_LEGs.csv"
data_LEGs = read.csv(source_LEGs, header = TRUE)

#Create dataframes for plot - consider mean vs median

md_total_1 <- mean(data_total[["Pos1_difference"]])
md_total_2 <- mean(data_total[["Pos2_difference"]])
md_total_3 <- mean(data_total[["Pos3_difference"]])
md_total_4 <- mean(data_total[["Pos4_difference"]])
md_total_5 <- mean(data_total[["Pos5_difference"]])

md_HEGs_1 <- mean(data_HEGs[["Pos1_difference"]])
md_HEGs_2 <- mean(data_HEGs[["Pos2_difference"]])
md_HEGs_3 <- mean(data_HEGs[["Pos3_difference"]])
md_HEGs_4 <- mean(data_HEGs[["Pos4_difference"]])
md_HEGs_5 <- mean(data_HEGs[["Pos5_difference"]])

md_LEGs_1 <- mean(data_LEGs[["Pos1_difference"]])
md_LEGs_2 <- mean(data_LEGs[["Pos2_difference"]])
md_LEGs_3 <- mean(data_LEGs[["Pos3_difference"]])
md_LEGs_4 <- mean(data_LEGs[["Pos4_difference"]])
md_LEGs_5 <- mean(data_LEGs[["Pos5_difference"]])

#Merge dataframes for plot

total <- c(md_total_1, md_total_2, md_total_3, md_total_4, md_total_5)
HEGs <- c(md_HEGs_1, md_HEGs_2, md_HEGs_3, md_HEGs_4, md_HEGs_5)
LEGs <- c(md_LEGs_1, md_LEGs_2, md_LEGs_3, md_LEGs_4, md_LEGs_5)

#Plots

par(mfrow=c(1,3), cex=1.05)

barplot(total, names.arg = c("+1","+2","+3","+4","+5"), xlab="Position", ylab="Median difference (ASC-absent - ASC-containing)", col = "#29B6F6", main="All genes", cex.main=2, cex.lab=1.5)

barplot(HEGs, names.arg = c("+1","+2","+3","+4","+5"), xlab="Position", ylab="Median difference (ASC-absent - ASC-containing)", col = "#29B6F6", main="HEGs", cex.main=2, cex.lab=1.5)

barplot(LEGs, names.arg = c("+1","+2","+3","+4","+5"), xlab="Position", ylab="Median difference (ASC-absent - ASC-containing)", col = "#29B6F6", main="LEGs", cex.main=2, cex.lab=1.5)

#Stats - wilcoxons
w1.total = wilcox.test(data_total[['Pos1_No_Stop']], data_total[['Pos1_Stop']], paired = TRUE, alternative = "greater")
w2.total = wilcox.test(data_total[['Pos2_No_Stop']], data_total[['Pos2_Stop']], paired = TRUE, alternative = "greater")
w3.total = wilcox.test(data_total[['Pos3_No_Stop']], data_total[['Pos3_Stop']], paired = TRUE, alternative = "greater")
w4.total = wilcox.test(data_total[['Pos4_No_Stop']], data_total[['Pos4_Stop']], paired = TRUE, alternative = "greater")
w5.total = wilcox.test(data_total[['Pos5_No_Stop']], data_total[['Pos5_Stop']], paired = TRUE, alternative = "greater")

w1.HEGs = wilcox.test(data_HEGs[['Pos1_No_Stop']], data_HEGs[['Pos1_Stop']], paired = TRUE, alternative = "greater")
w2.HEGs = wilcox.test(data_HEGs[['Pos2_No_Stop']], data_HEGs[['Pos2_Stop']], paired = TRUE, alternative = "greater")
w3.HEGs = wilcox.test(data_HEGs[['Pos3_No_Stop']], data_HEGs[['Pos3_Stop']], paired = TRUE, alternative = "greater")
w4.HEGs = wilcox.test(data_HEGs[['Pos4_No_Stop']], data_HEGs[['Pos4_Stop']], paired = TRUE, alternative = "greater")
w5.HEGs = wilcox.test(data_HEGs[['Pos5_No_Stop']], data_HEGs[['Pos5_Stop']], paired = TRUE, alternative = "greater")

w1.LEGs = wilcox.test(data_LEGs[['Pos1_No_Stop']], data_LEGs[['Pos1_Stop']], paired = TRUE, alternative = "greater")
w2.LEGs = wilcox.test(data_LEGs[['Pos2_No_Stop']], data_LEGs[['Pos2_Stop']], paired = TRUE, alternative = "greater")
w3.LEGs = wilcox.test(data_LEGs[['Pos3_No_Stop']], data_LEGs[['Pos3_Stop']], paired = TRUE, alternative = "greater")
w4.LEGs = wilcox.test(data_LEGs[['Pos4_No_Stop']], data_LEGs[['Pos4_Stop']], paired = TRUE, alternative = "greater")
w5.LEGs = wilcox.test(data_LEGs[['Pos5_No_Stop']], data_LEGs[['Pos5_Stop']], paired = TRUE, alternative = "greater")

#Stats - sign tests

library(BSDA)

s1.total = SIGN.test(data_total[['Pos1_No_Stop']], data_total[['Pos1_Stop']], paired = TRUE, alternative = "greater")
s2.total = SIGN.test(data_total[['Pos2_No_Stop']], data_total[['Pos2_Stop']], paired = TRUE, alternative = "greater")
s3.total = SIGN.test(data_total[['Pos3_No_Stop']], data_total[['Pos3_Stop']], paired = TRUE, alternative = "greater")
s4.total = SIGN.test(data_total[['Pos4_No_Stop']], data_total[['Pos4_Stop']], paired = TRUE, alternative = "greater")
s5.total = SIGN.test(data_total[['Pos5_No_Stop']], data_total[['Pos5_Stop']], paired = TRUE, alternative = "greater")

s1.HEGs = SIGN.test(data_HEGs[['Pos1_No_Stop']], data_HEGs[['Pos1_Stop']], paired = TRUE, alternative = "greater")
s2.HEGs = SIGN.test(data_HEGs[['Pos2_No_Stop']], data_HEGs[['Pos2_Stop']], paired = TRUE, alternative = "greater")
s3.HEGs = SIGN.test(data_HEGs[['Pos3_No_Stop']], data_HEGs[['Pos3_Stop']], paired = TRUE, alternative = "greater")
s4.HEGs = SIGN.test(data_HEGs[['Pos4_No_Stop']], data_HEGs[['Pos4_Stop']], paired = TRUE, alternative = "greater")
s5.HEGs = SIGN.test(data_HEGs[['Pos5_No_Stop']], data_HEGs[['Pos5_Stop']], paired = TRUE, alternative = "greater")

s1.LEGs = SIGN.test(data_LEGs[['Pos1_No_Stop']], data_LEGs[['Pos1_Stop']], paired = TRUE, alternative = "greater")
s2.LEGs = SIGN.test(data_LEGs[['Pos2_No_Stop']], data_LEGs[['Pos2_Stop']], paired = TRUE, alternative = "greater")
s3.LEGs = SIGN.test(data_LEGs[['Pos3_No_Stop']], data_LEGs[['Pos3_Stop']], paired = TRUE, alternative = "greater")
s4.LEGs = SIGN.test(data_LEGs[['Pos4_No_Stop']], data_LEGs[['Pos4_Stop']], paired = TRUE, alternative = "greater")
s5.LEGs = SIGN.test(data_LEGs[['Pos5_No_Stop']], data_LEGs[['Pos5_Stop']], paired = TRUE, alternative = "greater")
