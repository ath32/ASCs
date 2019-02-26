### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'
#source = '2_HEGs_fasta'
#source = '2_LEGs_fasta'

### For specific gene groups:
#source = '2_FASTA_Eubacteria_cds_TT11_TAA'
#source = '2_FASTA_Eubacteria_cds_TT11_TGA'
#source = '2_FASTA_Eubacteria_cds_TT11_TAG'

### FUNCTIONS ###

def main(source):

    total_csv = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            path = os.path.join(source, f)
            raw = open(path).read()

            #Obtain genome information
            a = raw.split('>')
            a1 = a[1].split(';')

            accession = a1[0]
            gc = float(a1[1].replace('GC', '').replace('=', ''))
            gc3 = float(a1[2].replace('GC3', '').replace('=', ''))

            if '>' in raw:

                genes = raw.strip().split('>')
                genes_c = list(filter(None, genes))

                #Get genes without an ASC at position 1
                no_asc_pos1 = get_ASC_absent(genes_c)

                #Get TNN frequency at position 1
                pos1_freq = get_TNN_pos1(no_asc_pos1)

                #Get average TNN frequency at positions 1 - 6
                utr_freq = get_TNN_utr(no_asc_pos1)

                #Create CSV line and append to total
                CSV_line = [accession, pos1_freq, utr_freq, gc, gc3]
                total_csv.append(CSV_line)

    headers = ["Accession", "F_p1", "F_utr", "GC", "GC3"]

    get_CSV(total_csv, headers)


def get_ASC_absent(chunks):

    #Create empty list to contain genes without an ASC at position +1
    list = []

    for i in chunks:

        #Get UTRs
        j = i.split('\n')
        utr = j[1]
        codons = [utr[i:i+3] for i in range(0, len(utr), 3)]

        #If codon 1 isn't a stop, add to list
        if codons[1] != 'taa' and codons[1] != 'tga' and codons[1] != 'tag':
            list.append(utr)

    return list


def get_TNN_pos1(utrs):

    #Set counter
    TNN_counter = 0

    for i in utrs:

        #Add 1 to counter if there is a 4th site T
        if i[3] == 't':
            TNN_counter += 1

    #Get frequency
    freq = TNN_counter / len(utrs)

    return freq


def get_TNN_utr(utrs):

    t_starting = 0
    n = 0

    #TNN codon list (minus stops)
    tnn_list = ['ttt', 'ttc', 'tta', 'ttg', 'tct', 'tcc', 'tca', 'tcg',
                'tat', 'tac', 'tgt', 'tgc', 'tgg']

    for j in utrs:

        #Get codons, ignoring the first one (primary stop)
        codons = [j[i:i+3] for i in range(0, len(j), 3)]
        utr_codons = codons[1:]

        for k in utr_codons:

            if k in tnn_list:
                t_starting += 1

        #Get total amount of codons tested
        n += len(utr_codons)


    freq = t_starting / n

    return freq


def get_CSV(total, headers):

    filename = "12.1_additional_1_v2.csv"
    #filename = "7.6_TAA.csv"
    #filename = "7.6_TGA.csv"
    #filename = "7.6_TAG.csv"
    subdir = "4_Outputs/CSVs"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in total:
            writer.writerow(j)


### RUN ###

if __name__ == '__main__':
    main(source)
