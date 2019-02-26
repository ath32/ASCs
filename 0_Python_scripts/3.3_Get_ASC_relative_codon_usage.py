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

                #Access each chunk
                genes = raw.strip().split('>')
                genes_c = list(filter(None, genes))

                #Only use genomes where there is more than 100 genes available
                #This is more appropriate when considering HEGs and LEGs
                if len(genes_c) >= 100:

                    #Obtain list of genes terminating with taa, tga, tag
                    taa_list, tga_list, tag_list = get_primary_lists(genes_c)

                    #Only consider genes that have at least one ASC at the primary stop
                    if check_genome(taa_list, tga_list, tag_list) == True:

                        #Calculate the number of ASCs for each list and what codon the ASC is
                        taa_taa_p, taa_tga_p, taa_tag_p = get_ASCs_proportions(taa_list)
                        tga_taa_p, tga_tga_p, tga_tag_p = get_ASCs_proportions(tga_list)
                        tag_taa_p, tag_tga_p, tag_tag_p = get_ASCs_proportions(tag_list)

                        csv_line = [accession, taa_taa_p, taa_tga_p, taa_tag_p, tga_taa_p, tga_tga_p, tga_tag_p,
                        tag_taa_p, tag_tga_p, tag_tag_p, GC, GC3]
                        csv_total.append(csv_line)

    headers = ['Accession', 'TAA-term TAA', 'TAA-term TGA', 'TAA-term TAG',
        'TGA-term TAA', 'TGA-term TGA', 'TGA-term TAG',
        'TAG-term TAA', 'TAG-term TGA', 'TAG-term TAG',
        'GC', 'GC3']

    create_csv(headers, csv_total)


def get_primary_lists(chunks):

    taa = []
    tga = []
    tag = []

    #Access UTR + primary stop
    for i in chunks:
        j = i.split('\n')
        utr = j[1]
        pri = utr[0:3]

    #Create list for each primary codon
        if pri == 'taa':
            taa.append(utr)
        elif pri == 'tga':
            tga.append(utr)
        elif pri == 'tag':
            tag.append(utr)

    return taa, tga, tag


def check_genome(taa, tga, tag):

    taa_asc = 0
    tga_asc = 0
    tag_asc = 0

    for i in taa:
        if i[3:6] == 'taa' or i[3:6] == 'tga' or i[3:6] == 'tag':
            taa_asc += 1
    for j in tga:
        if j[3:6] == 'taa' or j[3:6] == 'tga' or j[3:6] == 'tag':
            tga_asc += 1
    for k in tag:
        if k[3:6] == 'taa' or k[3:6] == 'tga' or k[3:6] == 'tag':
            tag_asc += 1

    if taa_asc > 0 and tga_asc > 0 and tag_asc > 0:
        outcome = True
    else:
        outcome = False

    return outcome


def get_ASCs_proportions(list):

    #Set counters for stops
    all_stops = 0
    taa = 0
    tga = 0
    tag = 0

    for utr in list:

        #Split into codons
        codons = [utr[i:i+3] for i in range(0, len(utr), 3)]

        #Calculate total ASCs in each UTR (only considering first position)
        if codons[1] == 'taa':
            all_stops += 1
            taa += 1
        elif codons[1] == 'tga':
            all_stops += 1
            tga += 1
        elif codons[1] == 'tag':
            all_stops += 1
            tag += 1

    taa_f = taa / all_stops
    tga_f = tga / all_stops
    tag_f = tag / all_stops

    return taa_f, tga_f, tag_f


def create_csv(headers, csv_total):

    ''' Output results to CSV file '''

    filename = "3.4_ASC_proportions_all.csv"
    #filename = "3.4_ASC_proportions_HEGs.csv"
    #filename = "3.4_ASC_proportions_LEGs.csv"
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
