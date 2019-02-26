### SOURCE

source  <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.9_all.csv"
data = read.csv(source, header = TRUE)

#Stats - compares TGA to other codons with A at the third site

w = wilcox.test(data[['TGAT_f']], data_all[['codonT_f']], paired = TRUE, alternative = 'g')


