### IMPORTS ###

import os
import numpy as np
import csv
import scipy.stats as st


### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'
#source = '2_HEGs_fasta'
#source = '2_LEGs_fasta'

def main(source):

    CSV_total = []

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

            #Access sequences and split into primary stop groups
            chunks = list(filter(None, a))
            taa_genes, tga_genes, tag_genes = get_t_lists(chunks)

            if len(taa_genes) > 0 and len(tga_genes) > 0 and len(tag_genes) > 0:

                #Get fifth site frequencies
                taa_a, taa_t, taa_g, taa_c = get_fifth_sites(taa_genes)
                tga_a, tga_t, tga_g, tga_c = get_fifth_sites(tga_genes)
                tag_a, tag_t, tag_g, tag_c = get_fifth_sites(tag_genes)

                #Create CSV line
                CSV_line = [accession, taa_a, taa_t, taa_g, taa_c, tga_a, tga_t, tga_g, tga_c, tag_a, tag_t, tag_g, tag_c, gc, gc3]

                CSV_total.append(CSV_line)

    headers = ['Accession', 'TAA_fifthA', 'TAA_fifthT', 'TAA_fifthG', 'TAA_fifthC',
        'TGA_fifthA', 'TGA_fifthT', 'TGA_fifthG', 'TGA_fifthC',
        'TAG_fifthA', 'TAG_fifthT', 'TAG_fifthG', 'TAG_fifthC',
        'GC', 'GC']

    get_CSV(CSV_total, headers)


def get_t_lists(chunks):

    taa = []
    tga = []
    tag = []

    for i in chunks:

        #Access utr
        b = i.split('\n')
        utr = b[1]

        if utr[3] == 't':

        #Split into groups
            if utr[0:3] == 'taa':
                taa.append(utr)
            elif utr[0:3] == 'tga':
                tga.append(utr)
            elif utr[0:3] == 'tag':
                tag.append(utr)

    return taa, tga, tag


def get_fifth_sites(utrs):

    total_genes = len(utrs)
    fifth_a = 0
    fifth_t = 0
    fifth_g = 0
    fifth_c = 0

    for i in utrs:

        #Get fourth site
        fifth = i[4]

        if fifth == 'a':
            fifth_a += 1
        if fifth == 't':
            fifth_t += 1
        if fifth == 'g':
            fifth_g += 1
        if fifth == 'c':
            fifth_c += 1

    freq_a = fifth_a / total_genes
    freq_t = fifth_t / total_genes
    freq_g = fifth_g / total_genes
    freq_c = fifth_c / total_genes

    return freq_a, freq_t, freq_g, freq_c


def get_CSV(total, headers):

    filename = "7.5_fifth_site_all.csv"
    #filename = "7.5_fifth_site_hegs.csv"
    #filename = "7.5_fifth_site_legs.csv"
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
