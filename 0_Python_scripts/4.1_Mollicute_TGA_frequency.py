### IMPORTS ###
import os
import numpy as np
import csv
import time

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_mollicutes'

### FUNCTIONS ###

def main(source):

    csv_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            path = os.path.join(source, f)
            raw = open(path).read()

            #Check file isn't empty
            if '>' in raw:

                #Obtain accession
                a = raw.split('>')
                a1 = a[1].split(';')
                accession = a1[0]

                #Check number of genes
                genes = raw.strip().split('>')
                genes_c = list(filter(None, genes))

                #Only use genomes where there is more than 100 genes available
                #This is more appropriate when considering HEGs and LEGs
                if len(genes_c) >= 100:

                    #Check translation transl_table
                    transl_table = a1[1]

                    ### UNHASH BELOW TO CHOOSE BETWEEN TT11 AND TT4 ###
                    #if transl_table == ' 11':
                    if transl_table == ' 4':

                        #Obtain GC and GC3
                        a2 = a1[2].strip()
                        a3 = a2.split("=")
                        GC = float(a3[1])

                        a4 = a1[3].strip()
                        a5 = a4.split("=")
                        GC3 = float(a5[1])

                        #Split into sections for use in functions
                        genes = raw.strip().split('>')
                        genes_c = list(filter(None, genes))

                        #Obtain UTR list and generate total UTR string, which may be useful later
                        utr_list, total_utr = get_utr_stuff(genes_c)

                        #Get observed stop frequencies
                        stop_counts, n_utrs, OF_array, OF_list = get_observed(utr_list)
                        OF_list.pop(0)

                        #CSV creation
                        nested = [[accession], OF_list, [GC], [GC3]]
                        output_flattened = [item for sublist in nested for item in sublist]

                        csv_total.append(output_flattened)


    #Set headers for the CSV
    headers = ["Accession", "P1_Observed", "P2_Observed", "P3_Observed", "P4_Observed", "P5_Observed", "P6_Observed", "Genomic_GC", "Genomic_GC3"]

    create_csv(headers, csv_total)


def get_utr_stuff(list):

    ''' Obtain list containing all utr sequences and string of total utr '''

    utr_list = []
    total_utr = ''

    #For each gene in the FASTA file...
    for i in list:

        #First generate UTR list, subnested for each codon
        a = i.split("\n")
        utr_seq = a[1]
        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        utr_list.append(codon_seq)

        #Now create total UTR list which will be useful later
        total_utr += utr_seq

    return utr_list, total_utr


def get_observed(list):

    ''' Calculate ASC frequencies at each position '''

    codons = [0, 0, 0, 0, 0, 0, 0]
    total_utrs = len(list)

    for i in list:

        if i[0] == 'tga': #or i[0] == 'tga' or i[0] == 'taa':
            codons[0] += 1
        if i[1] == 'tga': #or i[1] == 'tga' or i[1] == 'taa':
            codons[1] += 1
        if i[2] == 'tga': #or i[2] == 'tga' or i[2] == 'taa':
            codons[2] += 1
        if i[3] == 'tga': #or i[3] == 'tga' or i[3] == 'taa':
            codons[3] += 1
        if i[4] == 'tga': #or i[4] == 'tga' or i[4] == 'taa':
            codons[4] += 1
        if i[5] == 'tga': #or i[5] == 'tga' or i[5] == 'taa':
            codons[5] += 1
        if i[6] == 'tga': #or i[6] == 'tga' or i[6] == 'taa':
            codons[6] += 1

    frequencies_array = np.array(codons) / len(list)
    frequencies_list = frequencies_array.tolist()

    return codons, total_utrs, frequencies_array, frequencies_list

def create_csv(headers, csv_total):

    ''' Create output files '''

    filename = "4.2_TT4_mollicute_ASCs.csv"
    #filename = "4.2_TT11_mollicute_ASCs.csv"
    subdir = "4_Outputs/CSVs"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main(source)
