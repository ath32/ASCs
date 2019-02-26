''' Calculate ASC frequencies for all groups of genes '''

### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'
#source = '2_HEGs_fasta'#
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

                zero_taa, zero_tga, zero_tag = get_rf_zero(genes_c)
                one_taa, one_tga, one_tag = get_rf_one(genes_c)
                two_taa, two_tga, two_tag = get_rf_two(genes_c)

                csv_line = [accession, zero_taa, zero_tga, zero_tag, one_taa, one_tga, one_tag, two_taa, two_tga, two_tag, GC, GC3]
                csv_total.append(csv_line)

    headers = ['Accession', '0_TAA', '0_TGA', '0_TAG', '1_TAA', '1_TGA', '1_TAG', '2_TAA', '2_TGA', '2_TAG', 'GC', 'GC3']

    create_csv(headers, csv_total)



def get_rf_zero(chunks):

    n = 0
    taa = 0
    tga = 0
    tag = 0

    for i in chunks:
        b = i.split('\n')
        utr_seq = b[1]
        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        useful = codon_seq[1:7]

        for j in useful:

            if j == 'taa':
                n += 1
                taa += 1
            elif j == 'tga':
                n += 1
                tga += 1
            elif j == 'tag':
                n += 1
                tag += 1

    taa_f = taa / n
    tga_f = tga / n
    tag_f = tag / n

    return taa_f, tga_f, tag_f


def get_rf_one(chunks):

    n = 0
    taa = 0
    tga = 0
    tag = 0

    for i in chunks:
        b = i.split('\n')
        utr_seq = b[1]

        #Move along one (ignore primary stop and first base of utr)
        shift = utr_seq[4:]

        #Convert to codons
        codon_seq = [shift[i:i+3] for i in range(0, len(shift), 3)]

        #Only the first 5 codons are valid within our chosen IGR range
        useful = codon_seq[0:5]

        for j in useful:

            if j == 'taa':
                n += 1
                taa += 1
            elif j == 'tga':
                n += 1
                tga += 1
            elif j == 'tag':
                n += 1
                tag += 1

    taa_f = taa / n
    tga_f = tga / n
    tag_f = tag / n

    return taa_f, tga_f, tag_f


def get_rf_two(chunks):

    n = 0
    taa = 0
    tga = 0
    tag = 0

    for i in chunks:
        b = i.split('\n')
        utr_seq = b[1]

        #Move along one (ignore primary stop and first base of utr)
        shift = utr_seq[5:]

        #Convert to codons
        codon_seq = [shift[i:i+3] for i in range(0, len(shift), 3)]

        #Only the first 5 codons are valid within our chosen IGR range
        useful = codon_seq[0:5]

        for j in useful:

            if j == 'taa':
                n += 1
                taa += 1
            elif j == 'tga':
                n += 1
                tga += 1
            elif j == 'tag':
                n += 1
                tag += 1

    taa_f = taa / n
    tga_f = tga / n
    tag_f = tag / n

    return taa_f, tga_f, tag_f


def create_csv(headers, csv_total):

    ''' Output results to CSV file '''

    filename = "12.1_OSCs.csv"
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
