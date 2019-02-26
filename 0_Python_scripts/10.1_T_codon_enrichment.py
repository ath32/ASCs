### IMPORTS ###

import os
import csv
import numpy as np

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'
#source = '2_HEGs_fasta'
#source = '2_LEGs_fasta'
#source = '2_FASTA_Eubacteria_cds_TT11_TAA'
#source = '2_FASTA_Eubacteria_cds_TT11_TGA'
#source = '2_FASTA_Eubacteria_cds_TT11_TAG'

### FUNCTIONS ###

def main():

    CSV_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            path = os.path.join(source, f)
            raw = open(path).read()

            if '>' in raw:

                #Split into gene list for use in functions
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

                #Get list
                utr_list = get_utrs(genes_c)

                if len(utr_list) > 0:

                    #Get position 1 frequencies
                    nested_1, codons = get_frequencies(utr_list, 1)

                    #Get positions 3-6
                    nested_3, codons = get_frequencies(utr_list, 3)
                    nested_4, codons = get_frequencies(utr_list, 4)
                    nested_5, codons = get_frequencies(utr_list, 5)
                    nested_6, codons = get_frequencies(utr_list, 6)

                    #Calculate average frequencies for 3-6
                    averages = get_averages(nested_3, nested_4, nested_5, nested_6)

                    #Only use genomes where average > 0
                    if zero_check(averages) == True:

                        #Combine data in format ready for CSV, and calculate differences
                        freq_pos1 = [i[1] for i in nested_1]
                        combined = [freq_pos1, averages]
                        output = [list(i) for i in zip(*combined)]

                        for i in output:

                            ratio = (i[0] / i[1]) - 1
                            i.append(ratio)


                        #Flatten to create CSV line for each genome
                        output_flattened = [item for sublist in output for item in sublist]

                        #Add accession, GC and GC3
                        output_with_info = [[accession], output_flattened, [GC], [GC3]]
                        CSV_line = [item for sublist in output_with_info for item in sublist]

                        CSV_total.append(CSV_line)


    headers = ['Accession', 'taa_1', 'taa_3-6', 'taa_ratio',
            'tga_1', 'tga_3-6', 'tga_ratio',
            'tag_1', 'tag_3-6', 'tag_ratio',
            'ttt_1', 'ttt_3-6', 'ttt_ratio',
            'tta_1', 'tta_3-6', 'tta_ratio',
            'ttc_1', 'ttc_3-6', 'ttc_ratio',
            'ttg_1', 'ttg_3-6', 'ttg_ratio',
            'tat_1', 'tat_3-6', 'tat_ratio',
            'tac_1', 'tac_3-6', 'tac_ratio',
            'tca_1', 'tca_3-6', 'tca_ratio',
            'tct_1', 'tct_3-6', 'tct_ratio',
            'tcc_1', 'tcc_3-6', 'tcc_ratio',
            'tcg_1', 'tcg_3-6', 'tcg_ratio',
            'tgt_1', 'tgt_3-6', 'tgt_ratio',
            'tgc_1', 'tgc_3-6', 'tgc_ratio',
            'tgg_1', 'tgg_3-6', 'tgg_ratio',
            'GC', 'GC3']

    create_csv(headers, CSV_total)

def get_utrs(chunks):

    utr_list = []

    for i in chunks:

        #Obtain UTRs
        split_genes = i.split('\n')
        utr = split_genes[1]
        codon_seq = [utr[i:i+3] for i in range(0, len(utr), 3)]
        utr_list.append(codon_seq)

    return utr_list


def get_frequencies(utr_list, position):

    #List of T starting codons
    codons = ['taa', 'tga', 'tag', 'ttt', 'tta', 'ttc', 'ttg',
            'tat', 'tac', 'tca', 'tct', 'tcc', 'tcg', 'tgt', 'tgc',
            'tgg']

    #Nest list in preparation of codons
    nest = [[i] for i in codons]

    #For each codon option...
    for codon in nest:

        #Add counter for each codon
        codon.append(0)

        #For each UTR, add to the counter if the codon is present
        for utr in utr_list:
            if utr[position] == codon[0]:
                codon[1] += 1

    #Convert raw counts to frequency
    output = [[i, j/len(utr_list)] for [i, j] in nest]

    return output, codons


def get_averages(nested_3, nested_4, nested_5, nested_6):

    #Get lists of just frequencies
    three = [i[1] for i in nested_3]
    four = [i[1] for i in nested_4]
    five= [i[1] for i in nested_5]
    six = [i[1] for i in nested_6]

    raw = [three, four, five, six]

    #Re-nest by codon frequency
    nested = [list(i) for i in zip(*raw)]

    #Calculate means
    means = [np.mean(i) for i in nested]

    return means

def zero_check(averages):

    if all(i > 0 for i in averages):
        outcome = True
    else:
        outcome = False

    return outcome


def create_csv(headers, csv_total):

    filename = "10.4_T_codons_all.csv"
    #filename = "10.4_T_codons_hegs.csv"
    #filename = "10.4_T_codons_legs.csv"
    #filename = "10.4_T_codons_taa.csv"
    #filename = "10.4_T_codons_tga.csv"
    #filename = "10.4_T_codons_tag.csv"
    subdir = "4_Outputs/CSVs"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
