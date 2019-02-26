### IMPORTS ###

import os
import csv
import scipy.stats as st

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_CDS_expression'

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

                #Split into sections for use in functions
                genes = raw.strip().split('>')
                genes_c = list(filter(None, genes))

                #Get Spearman's ranks
                rho, p = get_spearman(genes_c)

                csv_line = [accession, rho, p, GC, GC3]
                csv_total.append(csv_line)

    headers = ['Accession', 'Rho', 'P', 'GC', 'GC3']

    create_csv(headers, csv_total)


def get_spearman(chunks):

    expressions = []
    lengths = []

    for i in chunks:

        b = i.split('\n')
        b1 = b[0].split(';')
        expression = float(b1[6])
        length = len(b[1])

        #Add data to list for spearman's test
        expressions.append(expression)
        lengths.append(length)

    spearmans = st.spearmanr(expressions, lengths)
    rho = spearmans[0]
    p = spearmans[1]

    return rho, p


def create_csv(headers, csv_total):

    ''' Create output files '''

    filename = "14.2_expression_length.csv"
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
