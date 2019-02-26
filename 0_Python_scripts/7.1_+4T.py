### IMPORTS ###

import os
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'
#source = '2_HEGs_fasta'
#source = '2_LEGs_fasta'

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

            #Access sequences and split into two groups - +4T or otherwise, then calculate frequency of +4T genes
            chunks = list(filter(None, a))
            present_list, absent_list, total = get_4T_lists(chunks)
            present_freq = len(present_list) / total

            #Get frequency of T starting codons in sequences
            t_starting_codon_f = get_t_codons(chunks)

            #Create CSV line and append to total
            CSV_line = [accession, present_freq, t_starting_codon_f, gc, gc3]
            total_csv.append(CSV_line)

    headers = ["Accession", "F_4T_genes", "F_Tstart_codons",
            "GC", "GC3"]

    get_CSV(total_csv, headers)


def get_4T_lists(chunks):

    present = []
    absent = []

    total = len(chunks)

    for i in chunks:
        b = i.split('\n')
        utr = b[1]

        if utr[3] == 't':
            present.append(utr)
        else:
            absent.append(utr)

    return present, absent, total


def get_t_codons(chunks):

    t_starting = 0
    n = 0

    for i in chunks:
        b = i.split('\n')
        utr = b[1]

        #Get codons (and exclude the primary stop codon)
        codons = [utr[i:i+3] for i in range(0, len(utr), 3)]
        utr_codons = codons[1:]

        #Get raw counts of T starting codons
        for j in utr_codons:

            if j[0] == 't':
                t_starting += 1
            else:
                t_starting = t_starting

        #Get total amount of codons tested
        n += len(utr_codons)

    genome_utr_t_freq = t_starting / n

    return genome_utr_t_freq


def get_CSV(total, headers):

    filename = "7.1_+4T.csv"
    #filename = "7.1_+4T_hegs.csv"
    #filename = "7.1_+4T_legs.csv"
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
