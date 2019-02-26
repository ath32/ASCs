### IMPORTS ###

import os
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_HEGs_fasta'
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

            total_genes = len(chunks)

            taa = 0
            tga = 0
            tag = 0

            for i in chunks:
                b = i.split('\n')
                utr = b[1]

                if utr[0:3] == 'taa':
                    taa += 1
                if utr[0:3] == 'tga':
                    tga += 1
                if utr[0:3] == 'tag':
                    tag += 1

            taa_p = taa / total_genes * 100
            tga_p = tga / total_genes * 100
            tag_p = tag / total_genes * 100

            CSV_line = [accession, taa_p, tga_p, tag_p, gc, gc3]

            total_csv.append(CSV_line)

    headers = ['Accession', 'TAA%', 'TGA%', 'TAG%', 'GC', 'GC3']

    filename = "8.1_primary_stop_percentages_HEGs.csv"
    #filename = "8.1_primary_stop_percentages_LEGs.csv"
    subdir = "4_Outputs/CSVs"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in total_csv:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main(source)
