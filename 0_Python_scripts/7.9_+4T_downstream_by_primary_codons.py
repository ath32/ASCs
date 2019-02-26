### IMPORTS ###

import os
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'
#source = '2_FASTA_Eubacteria_cds_TT11_TAA'
#source = '2_FASTA_Eubacteria_cds_TT11_TGA'
#source = '2_FASTA_Eubacteria_cds_TT11_TAG'

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

            #Get gene chunks
            genes = raw.strip().split('>')
            genes_c = list(filter(None, genes))

            #Filter genomes to make sure we can use them
            if filters(genes_c) == True:

                #Get +4T frequency
                fourth_t_f = get_fourth_t(genes_c)

                #Get 3' ASC +4T frequency
                taat_f = get_ASC_ts(genes_c, 'taa')
                tgat_f = get_ASC_ts(genes_c, 'tga')
                tagt_f = get_ASC_ts(genes_c, 'tag')

                #Get general codon +4T frequency
                codon_f = get_codon_ts(genes_c)

            csv_line = [accession, fourth_t_f, taat_f, tgat_f, tagt_f, codon_f, gc ,gc3]
            total_csv.append(csv_line)

    headers = ['Accession', '+4T_f', 'TAAT_f', 'TGAT_f', 'TAGT_f', 'codonT_f', 'GC', 'GC3']

    create_csv(headers, total_csv)


def get_fourth_t(chunks):

    fourth_t = 0

    for i in chunks:
        j = i.split('\n')
        utr = j[1]

        fourth_nt = utr[3]

        if fourth_nt == 't':
            fourth_t += 1

    frequency = fourth_t / len(chunks)

    return frequency


def filters(chunks):

    tga = 0
    taa = 0
    tag = 0

    #Get all UTRs
    for i in chunks:
        j = i.split('\n')
        utr = j[1]

        #Split utrs into codons - and isolate positions 1 to 6
        codons = [utr[i:i+3] for i in range(0, len(utr), 3)]
        useful = codons[1:7]

        for j, k in enumerate(useful):
            if k == 'tga':
                if j < 5:
                    tga += 1
            if k == 'taa':
                if j < 5:
                    taa += 1
            if k == 'tag':
                if j < 5:
                    tag += 1

    if tga > 0 and taa > 0 and tag > 0:
        outcome = True
    else:
        outcome = False

    return outcome


def get_ASC_ts(chunks, codon):

    t = 0
    n = 0

    for i in chunks:
        j = i.split('\n')
        utr = j[1]

        #Split utrs into codons - and isolate positions 1 to 6
        codons = [utr[i:i+3] for i in range(0, len(utr), 3)]
        useful = codons[1:7]

        #For each codon, determine if there is an ASC
        for i,j in enumerate(useful):

            if j == codon:
                #If there is an ASC, make sure there is another codon to interrogate 3'
                if i < 5:
                    n += 1

                    #If there is another codon, check whether the first base is a T
                    if useful[i+1][0] == 't':
                        t += 1

    freq = t/n

    return freq


def get_codon_ts(chunks):

    t = 0
    n = 0

    codons = ['tta', 'tca', 'taa', 'tga', 'cta', 'cca', 'caa', 'cga', 'ata', 'aca', 'aga', 'aaa',
    'gta', 'gca', 'gaa', 'gga']

    for i in chunks:
        j = i.split('\n')
        utr = j[1]

        #Split utrs into codons - and isolate positions 1 to 6
        codons = [utr[i:i+3] for i in range(0, len(utr), 3)]
        useful = codons[1:7]

        #For each codon, determine if there is an ASC
        for i,j in enumerate(useful):

            if j in codons:

                #Make sure theres another codon to interrogate
                if i < 5:
                    n += 1

                    #If there is another codon, check whether the first base is a T
                    if useful[i+1][0] == 't':
                        t += 1

    freq = t/n

    return freq


def create_csv(headers, csv_total):

    ''' Output results to CSV file '''

    filename = "7.9_all.csv"
    #filename = "7.8_TAA.csv"
    #filename = "7.8_TGA.csv"
    #filename = "7.8_TAG.csv"
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
