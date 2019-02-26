### IMPORTS ###

import os
import csv
import numpy as np

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_mollicutes'
tt11_source = '2_FASTA_Eubacteria_cds_TT11'

### FUNCTIONS ###

def main():

    CSV_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            path = os.path.join(source, f)
            raw = open(path).read()

            if '>' in raw:

                #Split into gene list for use in functions
                genes = raw.strip().split('>')
                genes_c = list(filter(None, genes))

                #Obtain accession
                a = raw.split('>')
                a1 = a[1].split(';')
                accession = a1[0]

                #Obtain GC and GC3
                a2 = a1[2].strip()
                a3 = a2.split("=")
                GC = float(a3[1])

                a4 = a1[3].strip()
                a5 = a4.split("=")
                GC3 = float(a5[1])

                #Only use genomes where there is more than 100 genes available
                #This is more appropriate when considering HEGs and LEGs
                if len(genes_c) >= 100:

                    #Check translation transl_table
                    transl_table = a1[1]

                    if transl_table == ' 4':

                        #Now obtain ASC frequencies at each position
                        utr_list = get_list(genes_c)
                        m_array, m_list = get_TGA_freqs(utr_list)

                        #Get matched genomes
                        p_array, p_list = get_matched_genomes(tt11_source, GC3)

                        #Create CSV
                        CSV_line = [[accession], m_list, p_list, [GC], [GC3]]
                        CSV_flattened = [item for sublist in CSV_line for item in sublist]
                        print (CSV_flattened)
                        CSV_total.append(CSV_flattened)

    headers = ['Accession', 'MF_0', 'MF_1', 'MF_2', 'MF_3', 'MF_4', 'MF_5', 'MF_6',
            'PF_0', 'PF_1', 'PF_2', 'PF_3', 'PF_4', 'PF_5', 'PF_6',
            'GC', 'GC3'
    ]

    create_csv(headers, CSV_total)


def get_list(genes_c):

    ''' Split UTR sequence into codons '''

    utr_list = []

    for i in genes_c:
        b = i.split('\n')
        utr_seq = b[1]

        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        utr_list.append(codon_seq)

    return utr_list


def get_TGA_freqs(utr_list):

    ''' Get stop frequencies at each position '''

    codons = [0, 0, 0, 0, 0, 0, 0]

    for i in utr_list:

        if i[0] == 'tga':
            codons[0] += 1
        if i[1] == 'tga':
            codons[1] += 1
        if i[2] == 'tga':
            codons[2] += 1
        if i[3] == 'tga':
            codons[3] += 1
        if i[4] == 'tga':
            codons[4] += 1
        if i[5] == 'tga':
            codons[5] += 1
        if i[6] == 'tga':
            codons[6] += 1

    frequencies_array = np.array(codons) / len(utr_list)
    frequencies_list = frequencies_array.tolist()

    return frequencies_array, frequencies_list


def get_matched_genomes(tt11_source, GC3):

    #Create empty list, wil contain lists of ASC frequencies from each GC matched genome
    freq_list = []

    for root, dirs, filenames in os.walk(tt11_source):
        for f in filenames:

            path1 = os.path.join(tt11_source, f)
            raw = open(path1).read()

            #Split into gene list for use in functions
            genes = raw.strip().split('>')
            genes_c = list(filter(None, genes))

            #Obtain accession
            a = raw.split('>')
            a1 = a[1].split(';')
            tt11_accession = a1[0]

            #Obtain GC and GC3
            a2 = a1[1].strip()
            a3 = a2.split("=")
            tt_11GC = float(a3[1])

            a4 = a1[2].strip()
            a5 = a4.split("=")
            tt11_GC3 = float(a5[1])

            if (GC3 - 3.5) < tt11_GC3 < (GC3 + 3.5):

                utr_list = get_list(genes_c)
                genome_array, genome_list = get_TGA_freqs(utr_list)
                freq_list.append(genome_list)

    freq_array = np.array(freq_list)
    freq_mean = np.mean(freq_array, axis=0)
    to_list = freq_mean.tolist()

    return freq_array, to_list


def create_csv(headers, csv_total):

    filename = "4.3_GC_matched.csv"
    subdir = "4_Outputs/CSVs"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
