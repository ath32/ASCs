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

            #Get gene chunks
            genes = raw.strip().split('>')
            genes_c = list(filter(None, genes))

            #Get 3' ASC +4T frequency
            ASC_f = get_ASC_ts(genes_c)

            #Get 3' codon +4T frequency
            codon_f = get_codon_ts(genes_c)

            csv_line = [accession, ASC_f, codon_f, gc, gc3]
            total_csv.append(csv_line)

    headers = ['Accession', 'ASC_f', 'Codon_f', 'GC', 'GC3']

    create_csv(headers, total_csv)


def get_ASC_ts(chunks):

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

            if j == 'tag' or j == 'tga' or j == 'taa':
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

    for i in chunks:
        j = i.split('\n')
        utr = j[1]

        #Split utrs into codons - and isolate positions 1 to 6
        codons = [utr[i:i+3] for i in range(0, len(utr), 3)]
        useful = codons[1:7]

        #For each codon, determine if there is an ASC
        for i,j in enumerate(useful):

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

    filename = "13.1_downstream_+4t.csv"
    #filename = "13.1_downstream_+4t_hegs.csv"
    #filename = "13.1_downstream_+4t_legs.csv"
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
