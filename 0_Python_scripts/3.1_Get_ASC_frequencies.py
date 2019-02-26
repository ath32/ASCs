### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'
#source = '2_HEGs_fasta'
#source = '2_LEGs_fasta'

### FUNCTIONS ###

def main(source):

    csv_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            path = os.path.join(source, f)
            raw = open(path).read()

            if '>' in raw:

                genes = raw.strip().split('>')
                genes_c = list(filter(None, genes))

                #Only use genomes where there is more than 100 genes available
                #This is more appropriate when considering HEGs and LEGs
                if len(genes_c) >= 100:

                    #Obtain accession
                    a = raw.split('>')
                    a1 = a[1].split(';')
                    accession = a1[0]

                    #Obtain GC and GC3
                    a2 = a1[1].strip()
                    a3 = a2.split("=")
                    GC = float(a3[1])

                    a4 = a1[2].strip()
                    a5 = a4.split("=")
                    GC3 = float(a5[1])

                    #Split into sections for use in functions
                    genes = raw.strip().split('>')
                    genes_c = list(filter(None, genes))

                    #Split UTRs into codons
                    glist = get_list(genes_c)

                    #Calculate frequencies for each group of genes
                    stop_freq_a, stop_freq_l = get_stop_freq(glist)
                    taa_freq_a, taa_freq_l = get_freqs(glist, 'taa')
                    tga_freq_a, tga_freq_l = get_freqs(glist, 'tga')
                    tag_freq_a, tag_freq_l = get_freqs(glist, 'tag')

                    #Manipulate data for suitability in CSV file
                    raw_list = [stop_freq_l, taa_freq_l, tag_freq_l, tga_freq_l]
                    nested = [list(i) for i in zip(*raw_list)]

                    #Flatten list
                    output_flattened = [item for sublist in nested for item in sublist]

                    #Create master list, which will represents a row of the output CSV
                    csv_nested = [[accession], output_flattened, [GC], [GC3]]
                    csv_list = [item for sublist in csv_nested for item in sublist]
                    csv_total.append(csv_list)

                else:
                    continue

    #Set headers for the CSV
    headers = ["Accession", "P0_stop_f", "P0_TAA_freq", "P0_TAG_f", "P0_TGA_f",
        "P1_stop_f", "P1_TAA_freq", "P1_TAG_f", "P1_TGA_f",
        "P2_stop_f", "P2_TAA_freq", "P2_TAG_f", "P2_TGA_f",
        "P3_stop_f", "P3_TAA_freq", "P3_TAG_f", "P3_TGA_f",
        "P4_stop_f", "P4_TAA_freq", "P4_TAG_f", "P4_TGA_f",
        "P5_stop_f", "P5_TAA_freq", "P5_TAG_f", "P5_TGA_f",
        "P6_stop_f", "P6_TAA_freq", "P6_TAG_f", "P6_TGA_f",
        "Genomic_GC", "Genomic_GC3"]

    create_csv(headers, csv_total)


def get_list(genes_c):

    ''' Split UTR sequence into codons '''

    utr_list = []

    for i in genes_c:
        b = i.split('\n')
        utr_seq = b[1]

        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        utr_list.append(codon_seq)

    return utr_list


def get_stop_freq(utr_list):

    ''' Get stop frequencies at each position '''

    codons = [0, 0, 0, 0, 0, 0, 0]

    for i in utr_list:

        if i[0] == 'tag' or i[0] == 'tga' or i[0] == 'taa':
            codons[0] += 1
        if i[1] == 'tag' or i[1] == 'tga' or i[1] == 'taa':
            codons[1] += 1
        if i[2] == 'tag' or i[2] == 'tga' or i[2] == 'taa':
            codons[2] += 1
        if i[3] == 'tag' or i[3] == 'tga' or i[3] == 'taa':
            codons[3] += 1
        if i[4] == 'tag' or i[4] == 'tga' or i[4] == 'taa':
            codons[4] += 1
        if i[5] == 'tag' or i[5] == 'tga' or i[5] == 'taa':
            codons[5] += 1
        if i[6] == 'tag' or i[6] == 'tga' or i[6] == 'taa':
            codons[6] += 1

    frequencies_array = np.array(codons) / len(utr_list)
    frequencies_list = frequencies_array.tolist()

    return frequencies_array, frequencies_list


def get_freqs(list, codon):

    ''' Get frequencies for a specific codon '''

    codons = [0, 0, 0, 0, 0, 0, 0]

    for i in list:

        if i[0] == codon:
            codons[0] += 1
        if i[1] == codon:
            codons[1] += 1
        if i[2] == codon:
            codons[2] += 1
        if i[3] == codon:
            codons[3] += 1
        if i[4] == codon:
            codons[4] += 1
        if i[5] == codon:
            codons[5] += 1
        if i[6] == codon:
            codons[6] += 1

    frequencies_array = np.array(codons) / len(list)
    frequencies_list = frequencies_array.tolist()

    return frequencies_array, frequencies_list


def create_csv(headers, csv_total):

    ''' Output results to CSV file '''

    filename = "3.1_Frequencies_all.csv"
    #filename = "3.1_Frequencies_HEGs.csv"
    #filename = "3.1_Frequencies_LEGs.csv"
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
