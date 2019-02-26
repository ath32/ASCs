### SOURCE ### Unhash source file of interest

source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/13.1_downstream_+4t.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/13.1_downstream_+4t_hegs.csv"
#source <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/13.1_downstream_+4t_legs.csv"
data = read.csv(source, header = TRUE)

#Calculate differences
data[['diffs']] <- data[['ASC_f']] - data[['Codon_f']]

#Plot differences
par(mfrow=c(1,3))

plot(data[['GC3']], data[['ASC_f']], xlim = c(0, 100), xlab='GC3', ylab='+4T frequency', main='ASCs', pch=16, col='black')
plot(data[['GC3']], data[['Codon_f']], xlim = c(0, 100), xlab='GC3', ylab='+4T frequency', main='All codons', pch=16, col='black')
plot(data[['diffs']])

#Stats - compare ASCs to general non-stop codons
w = wilcox.test(data[['ASC_f']], data[['Codon_f']], paired = TRUE, alternative = 'l')

### Stats part 2 - compare ASCs to the primary stop
### Requires additional source - unhash the one that matches the source at the top of the script

source2 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.1_+4T.csv"
#source2 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.1_+4T_hegs.csv"
#source2 <- "/Users/atho/Documents/Project_ASCs/4_Outputs/CSVs/7.1_+4T_legs.csv"
data2 = read.csv(source2, header = TRUE)

ASCs <- data[['ASC_f']]
Primary <- data2[['F_4T_genes']]

w2 = wilcox.test(ASCs, Primary, paired = TRUE, alternative = 'l')

