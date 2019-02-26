### SOURCE ###

#Requires source for all genes (below)
source_all  <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.8_all.csv"
data_all = read.csv(source_all, header = TRUE)

#Also requires source for TAA, TGA or TAG-terminating genes, depending on the test you want to do. (unhash the source file of interest)
source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.8_TAA.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.8_TGA.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.8_TAG.csv"
data = read.csv(source, header = TRUE)

#Test 1 - compare stop codon variant as ASC to primary stop

w1 = wilcox.test(data_all[['TAAT_f']], data[['X.4T_f']], paired = TRUE, alternative = 'l')
#w1 = wilcox.test(data[['TGAT_f']], data[['X.4T_f']], paired = TRUE, alternative = 'l')
#w1 = wilcox.test(data[['TAGT_f']], data[['X.4T_f']], paired = TRUE, alternative = 'l')

#Test 2 - compare ASC variant to other non-stop codons

w2 = wilcox.test(data_all[['TAAT_f']], data_all[['codonT_f']], paired = TRUE, alternative = 'l')
#w2 = wilcox.test(data[['TGAT_f']], data_all[['codonT_f']], paired = TRUE, alternative = 'g')
#w2 = wilcox.test(data[['TAGT_f']], data_all[['codonT_f']], paired = TRUE, alternative = 'l')


#Test 3 - compare ASC variants to each other

w3 = wilcox.test(data_all[['TGAT_f']], data_all[['TAAT_f']], paired = TRUE, alternative = 'g')
w4 = wilcox.test(data_all[['TGAT_f']], data_all[['TAGT_f']], paired = TRUE, alternative = 'g')
w5 = wilcox.test(data_all[['TAAT_f']], data_all[['TAGT_f']], paired = TRUE, alternative = 'g')

### Analysis of expression groups 

source_hegs  <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.8_hegs.csv"
data_hegs = read.csv(source_hegs, header = TRUE)

source_legs  <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.8_legs.csv"
data_legs = read.csv(source_legs, header = TRUE)

w5 = wilcox.test(data_hegs[['TGAT_f']], data_legs[['TGAT_f']], paired = TRUE, alternative = 'two.sided')
w6 = wilcox.test(data_hegs[['TAAT_f']], data_legs[['TAAT_f']], paired = TRUE, alternative = 'g')
w7 = wilcox.test(data_hegs[['TAGT_f']], data_legs[['TAGT_f']], paired = TRUE, alternative = 'g')

a = mean(data_hegs[['TGAT_f']]) - mean(data_legs[['TGAT_f']])
b = mean(data_hegs[['TAAT_f']]) - mean(data_legs[['TAAT_f']])
c = mean(data_hegs[['TAGT_f']]) - mean(data_legs[['TAGT_f']])

d = median(data_hegs[['TGAT_f']]) - median(data_legs[['TGAT_f']])
e = median(data_hegs[['TAAT_f']]) - median(data_legs[['TAAT_f']])
f = median(data_hegs[['TAGT_f']]) - median(data_legs[['TAGT_f']])

