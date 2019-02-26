### IMPORTS ###

import os
import numpy as np
import csv
import time
import scipy.stats as st

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_HEGs_fasta'
#source = '2_LEGs_fasta'

### FUNCTIONS ###

def main(source):

    csv_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            path = os.path.join(source, f)
            raw = open(path).read()

            #Check file isn't empty
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

                #Split into sections for use in functions
                genes = raw.strip().split('>')
                genes_c = list(filter(None, genes))

                #Only consider genomes with over 100 candidate genes
                if len(genes_c) >= 100:

                    #Get lists of genes with each primary stop
                    taa_list, tga_list, tag_list, non_taa_list = get_primary_stop_lists(genes_c)

                    if len(taa_list) > 0 and  len(tga_list) > 0 and len(tag_list) > 0:

                        #Get stop frequencies
                        taa_array, taa_list = get_stop_freq(taa_list)
                        tga_array, tga_list = get_stop_freq(tga_list)
                        tag_array, tag_list = get_stop_freq(tag_list)
                        non_taa_array, non_taa_list = get_stop_freq(non_taa_list)


                        raw_list = [taa_list, non_taa_list, tag_list, tga_list]
                        nested = [list(i) for i in zip(*raw_list)]

                        #Flatten list
                        output_flattened = [item for sublist in nested for item in sublist]

                        #Create master list, which will represents a row of the output CSV
                        csv_nested = [[accession], output_flattened, [GC], [GC3]]
                        csv_list = [item for sublist in csv_nested for item in sublist]

                        csv_total.append(csv_list)

    #Set headers for the CSV
    headers = ["Accession", "P0_Primary_TAA", "P0_Primary_nonTAA", "P0_Primary_TAG", "P0_Primary_TGA",
        "P1_Primary_TAA", "P1_Primary_nonTAA", "P1_Primary_TAG", "P1_Primary_TGA",
        "P2_Primary_TAA", "P2_Primary_nonTAA", "P2_Primary_TAG", "P2_Primary_TGA",
        "P3_Primary_TAA", "P3_Primary_nonTAA", "P3_Primary_TAG", "P3_Primary_TGA",
        "P4_Primary_TAA", "P4_Primary_nonTAA", "P4_Primary_TAG", "P4_Primary_TGA",
        "P5_Primary_TAA", "P5_Primary_nonTAA", "P5_Primary_TAG", "P5_Primary_TGA",
        "P6_Primary_TAA", "P6_Primary_nonTAA", "P6_Primary_TAG", "P6_Primary_TGA",
        "Genomic_GC", "Genomic_GC3"]

    create_csv(headers, csv_total)


def get_primary_stop_lists(chunks):

    no_fourth_t = []

    taa = []
    tga = []
    tag = []
    non_taa = []

    for i in chunks:
        b = i.split('\n')
        utr_seq = b[1]

        if utr_seq[3] != 't':
            no_fourth_t.append(utr_seq)

    for j in no_fourth_t:
        codon_seq = [j[i:i+3] for i in range(0, len(j), 3)]

        if codon_seq[0] == 'taa':
            taa.append(codon_seq)
        elif codon_seq[0] == 'tga':
            tga.append(codon_seq)
            non_taa.append(codon_seq)
        elif codon_seq[0] == 'tag':
            tag.append(codon_seq)
            non_taa.append(codon_seq)

    return taa, tga, tag, non_taa


def get_stop_freq(list):

    codons = [0, 0, 0, 0, 0, 0, 0]

    for i in list:

        if i[0] == 'tag' or i[0] == 'tga' or i[0] == 'taa':
            codons[0] += 1
        if i[1] == 'tag' or i[1] == 'tga' or i[1] == 'taa':
            codons[1] += 1
        if i[2] == 'tag' or i[2] == 'tga' or i[2] == 'taa':
            codons[2] += 1
        if i[3] == 'tag' or i[3] == 'tga' or i[3] == 'taa':
            codons[3] += 1
        if i[4] == 'tag' or i[4] == 'tga' or i[4] == 'taa':
            codons[4] += 1
        if i[5] == 'tag' or i[5] == 'tga' or i[5] == 'taa':
            codons[5] += 1
        if i[6] == 'tag' or i[6] == 'tga' or i[6] == 'taa':
            codons[6] += 1

    frequencies_array = np.array(codons) / len(list)
    frequencies_array[np.isnan(frequencies_array)] = 0
    frequencies_list = frequencies_array.tolist()

    return frequencies_array, frequencies_list


def create_csv(headers, csv_total):

    filename = "5.2_Subset_diffs_HEGs.csv"
    #filename = "5.2_Subset_diffs_LEGs.csv"
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
