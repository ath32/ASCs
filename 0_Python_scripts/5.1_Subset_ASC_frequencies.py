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

                #Create new strings dependent on primary stop codon
                taa_genes = ''
                tag_genes = ''
                tga_genes = ''
                non_taa_genes = ''

                if len(genes_c) >= 100:

                    for i in genes_c:

                        #Obtain primary stop codons
                        c = i.split('\n')
                        c1 = c[1]
                        primary = c1[:3]

                        #Append to one of the three lists
                        if primary == 'taa':
                            taa_genes += ">" + i
                        elif primary == 'tag':
                            tag_genes += ">" + i
                            non_taa_genes += ">" + i
                        elif primary == 'tga':
                            tga_genes += ">" + i
                            non_taa_genes += ">" + i
                        else:
                            print ('no primary stop')

                else:
                    continue

                #Create lists for use in functions
                list_taa = taa_genes.split('>')
                list_tag = tag_genes.split('>')
                list_tga = tga_genes.split('>')
                list_non_taa = non_taa_genes.split('>')

                #Filter to remove empty entries and split into codons for functions
                taa_c = list(filter(None, list_taa))
                taa = get_list(taa_c)

                tag_c = list(filter(None, list_tag))
                tag = get_list(tag_c)

                tga_c = list(filter(None, list_tga))
                tga = get_list(tga_c)

                non_taa_c = list(filter(None, list_non_taa))
                non_taa = get_list(non_taa_c)

                #Calculate stop frequencies
                taa_freq_a, taa_freq_l = get_stop_freq(taa)
                tag_freq_a, tag_freq_l = get_stop_freq(tag)
                tga_freq_a, tga_freq_l = get_stop_freq(tga)
                non_taa_freq_a, non_taa_freq_l = get_stop_freq(non_taa)

                raw_list = [taa_freq_l, non_taa_freq_l, tag_freq_l, tga_freq_l]
                nested = [list(i) for i in zip(*raw_list)]

                #Flatten list
                output_flattened = [item for sublist in nested for item in sublist]

                #Create master list, which will represents a row of the output CSV
                csv_nested = [[accession], output_flattened, [GC], [GC3]]
                csv_list = [item for sublist in csv_nested for item in sublist]

                csv_total.append(csv_list)

            else:
                continue

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


def get_list(genes_c):

    utr_list = []

    for i in genes_c:
        b = i.split('\n')
        utr_seq = b[1]

        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        utr_list.append(codon_seq)

    return utr_list


def get_stop_freq(utr_list):

    codons = [0, 0, 0, 0, 0, 0, 0]

    for i in utr_list:

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

    frequencies_array = np.array(codons) / len(utr_list)
    frequencies_array[np.isnan(frequencies_array)] = 0
    frequencies_list = frequencies_array.tolist()

    return frequencies_array, frequencies_list

def get_freqs(list, codon):

    codons = [0, 0, 0, 0, 0, 0, 0]

    for i in list:

        if i[0] == codon:
            codons[0] += 1
        if i[1] == codon:
            codons[1] += 1
        if i[2] == codon:
            codons[2] += 1
        if i[3] == codon:
            codons[3] += 1
        if i[4] == codon:
            codons[4] += 1
        if i[5] == codon:
            codons[5] += 1
        if i[6] == codon:
            codons[6] += 1

    frequencies_array = np.array(codons) / len(list)
    frequencies_list = frequencies_array.tolist()

    return frequencies_array, frequencies_list

def create_csv(headers, csv_total):

    filename = "5.1_Subset_diffs_HEGs.csv"
    #filename = "5.1_Subset_diffs_LEGs.csv"
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
