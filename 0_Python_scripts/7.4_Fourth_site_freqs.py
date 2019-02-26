### IMPORTS ###

import os
import csv

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
            taa_genes, tga_genes, tag_genes = get_lists(chunks)

            if len(taa_genes) > 0 and len(tga_genes) > 0 and len(tag_genes) > 0:

                #Get fourth site frequencies
                taa_a, taa_t, taa_g, taa_c = get_fourth_sites(taa_genes)
                tga_a, tga_t, tga_g, tga_c = get_fourth_sites(tga_genes)
                tag_a, tag_t, tag_g, tag_c = get_fourth_sites(tag_genes)

                #Create CSV line
                CSV_line = [accession, taa_a, taa_t, taa_g, taa_c, tga_a, tga_t, tga_g, tga_c, tag_a, tag_t, tag_g, tag_c, gc, gc3]

                CSV_total.append(CSV_line)

    headers = ['Accession', 'TAA_FourthA', 'TAA_FourthT', 'TAA_FourthG', 'TAA_FourthC',
        'TGA_FourthA', 'TGA_FourthT', 'TGA_FourthG', 'TGA_FourthC',
        'TAG_FourthA', 'TAG_FourthT', 'TAG_FourthG', 'TAG_FourthC',
        'GC', 'GC']

    get_CSV(CSV_total, headers)


def get_lists(chunks):

    taa = []
    tga = []
    tag = []

    for i in chunks:

        #Access utr
        b = i.split('\n')
        utr = b[1]

        #Split into groups
        if utr[0:3] == 'taa':
            taa.append(utr)
        elif utr[0:3] == 'tga':
            tga.append(utr)
        elif utr[0:3] == 'tag':
            tag.append(utr)

    return taa, tga, tag


def get_fourth_sites(utrs):

    total_genes = len(utrs)
    fourth_a = 0
    fourth_t = 0
    fourth_g = 0
    fourth_c = 0

    for i in utrs:

        #Get fourth site
        fourth = i[3]

        if fourth == 'a':
            fourth_a += 1
        if fourth == 't':
            fourth_t += 1
        if fourth == 'g':
            fourth_g += 1
        if fourth == 'c':
            fourth_c += 1

    freq_a = fourth_a / total_genes
    freq_t = fourth_t / total_genes
    freq_g = fourth_g / total_genes
    freq_c = fourth_c / total_genes

    return freq_a, freq_t, freq_g, freq_c


def get_CSV(total, headers):

    filename = "7.4_fourth_site_all.csv"
    #filename = "7.4_fourth_site_hegs.csv"
    #filename = "7.4_fourth_site_legs.csv"
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
