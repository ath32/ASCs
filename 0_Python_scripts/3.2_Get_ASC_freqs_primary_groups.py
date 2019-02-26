''' Calculate ASC frequencies for all groups of genes '''

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

            #Ignore any empty files
            if '>' in raw:

                #Access each chunk
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

                    #Split intro three UTR groups - TAA, TGA, TAG-terminating genes
                    taa, tga, tag = get_utr_lists(genes_c)

                    #Get stop codon frequencies at each position
                    if len(taa) > 0 and len(tga) > 0 and len(tag) > 0:
                        taa_freqs, tga_freqs, tag_freqs = get_freqs(taa, tga, tag)

                        #Write output line for genome
                        output_line = [[accession], taa_freqs, tga_freqs, tag_freqs, [GC], [GC3]]
                        CSV_line = [item for sublist in output_line for item in sublist]
                        csv_total.append(CSV_line)

    headers = ["Accession", "TAA_0", "TAA_1", "TAA_2", "TAA_3", "TAA_4", "TAA_5", "TAA_6",
        "TGA_0", "TGA_1", "TGA_2", "TGA_3", "TGA_4", "TGA_5", "TGA_6",
        "TAG_0", "TAG_1", "TAG_2", "TAG_3", "TAG_4", "TAG_5", "TAG_6",
        "GC", "GC3"
    ]

    create_csv(headers, csv_total)


def get_utr_lists(chunks):

    #Set empty lists
    taa = []
    tga = []
    tag = []

    #First access the UTR sequence
    for i in chunks:
        j = i.split('\n')
        utr = j[1]

        #Next, split utrs into codons
        codons = [utr[i:i+3] for i in range(0, len(utr), 3)]

        #Categorize into three lists
        if codons[0] == 'taa':
            taa.append(codons)
        elif codons[0] == 'tga':
            tga.append(codons)
        elif codons[0] == 'tag':
            tag.append(codons)

    return taa, tga, tag


def get_freqs(taa, tga, tag):

    #Set stop counters
    taa_count = [0, 0, 0, 0, 0, 0, 0]
    tga_count = [0, 0, 0, 0, 0, 0, 0]
    tag_count = [0, 0, 0, 0, 0, 0, 0]

    #TAA
    for i in taa:
        if i[0] == 'taa' or i[0] == 'tga' or i[0] == 'tag':
            taa_count[0] += 1
        if i[1] == 'taa' or i[1] == 'tga' or i[1] == 'tag':
            taa_count[1] += 1
        if i[2] == 'taa' or i[2] == 'tga' or i[2] == 'tag':
            taa_count[2] += 1
        if i[3] == 'taa' or i[3] == 'tga' or i[3] == 'tag':
            taa_count[3] += 1
        if i[4] == 'taa' or i[4] == 'tga' or i[4] == 'tag':
            taa_count[4] += 1
        if i[5] == 'taa' or i[5] == 'tga' or i[5] == 'tag':
            taa_count[5] += 1
        if i[6] == 'taa' or i[6] == 'tga' or i[6] == 'tag':
            taa_count[6] += 1

    #TGA
    for i in tga:
        if i[0] == 'taa' or i[0] == 'tga' or i[0] == 'tag':
            tga_count[0] += 1
        if i[1] == 'taa' or i[1] == 'tga' or i[1] == 'tag':
            tga_count[1] += 1
        if i[2] == 'taa' or i[2] == 'tga' or i[2] == 'tag':
            tga_count[2] += 1
        if i[3] == 'taa' or i[3] == 'tga' or i[3] == 'tag':
            tga_count[3] += 1
        if i[4] == 'taa' or i[4] == 'tga' or i[4] == 'tag':
            tga_count[4] += 1
        if i[5] == 'taa' or i[5] == 'tga' or i[5] == 'tag':
            tga_count[5] += 1
        if i[6] == 'taa' or i[6] == 'tga' or i[6] == 'tag':
            tga_count[6] += 1

    #TAG
    for i in tag:
        if i[0] == 'taa' or i[0] == 'tga' or i[0] == 'tag':
            tag_count[0] += 1
        if i[1] == 'taa' or i[1] == 'tga' or i[1] == 'tag':
            tag_count[1] += 1
        if i[2] == 'taa' or i[2] == 'tga' or i[2] == 'tag':
            tag_count[2] += 1
        if i[3] == 'taa' or i[3] == 'tga' or i[3] == 'tag':
            tag_count[3] += 1
        if i[4] == 'taa' or i[4] == 'tga' or i[4] == 'tag':
            tag_count[4] += 1
        if i[5] == 'taa' or i[5] == 'tga' or i[5] == 'tag':
            tag_count[5] += 1
        if i[6] == 'taa' or i[6] == 'tga' or i[6] == 'tag':
            tag_count[6] += 1


    #Get length of gene lists
    taa_l = len(taa)
    tga_l = len(tga)
    tag_l = len(tag)

    #Calculate frequencies
    taa_f = np.array(taa_count) / taa_l
    tga_f = np.array(tga_count) / tga_l
    tag_f = np.array(tag_count) / tag_l

    #Convert back to lists
    taa_fl = taa_f.tolist()
    tga_fl = tga_f.tolist()
    tag_fl = tag_f.tolist()

    return taa_fl, tga_fl, tag_fl


def create_csv(headers, csv_total):

    ''' Output results to CSV file '''

    filename = "3.3_Freq_pri_all.csv"
    #filename = "3.3_Freq_pri_HEGs.csv"
    #filename = "3.3_Freq_pri_LEGs.csv"
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
