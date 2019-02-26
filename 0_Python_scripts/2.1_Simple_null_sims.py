### IMPORTS ###

import os
import numpy as np
import csv
import time
import scipy.stats as st

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '1_FASTA_Eubacteria_cds_TT11'

### FUNCTIONS ###

def folder_parse(source):

    csv_total = []

    for root, dirs, filenames in os.walk(source):

        for f in filenames:
            path = os.path.join(source, f)
            raw = open(path).read()

            #Start timer
            start_time = time.time()

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

                #Obtain UTR list and generate total UTR string, which is useful for later functions
                utr_list, total_utr = get_utr_stuff(genes_c)

                #Obtain observed frequencies of stop codons at each codon position
                stop_counts, n_utrs, OF_array, OF_list = get_freq(utr_list)

                #Array in case its needed
                stop_counts_arr = np.array(stop_counts)

                #Simulate 10,000 UTRs according to dinucleotide content to generate expected stop codon frequencies
                EF_array, EF_list = get_simulations(utr_list, total_utr, start_time)

                #Create nested list, containing all elements needed for the binomial test
                stats_list = [stop_counts, [n_utrs, n_utrs, n_utrs, n_utrs, n_utrs, n_utrs, n_utrs], EF_list]
                nested = [list(i) for i in zip(*stats_list)]

                #Calculate and add observed frequency, standard deviation, z-score and binomial p-value to the nested list
                for i in nested:
                    observed_freq = i[0] / i[1]
                    i.append(observed_freq)

                    standard_dev = (10000 * (i[2]) * (1-i[2])) ** 0.5
                    i.append(standard_dev)

                    z_score = (i[0] - (i[1] * i[2])) / i[4]
                    i.append(z_score)

                    p = st.binom_test(i[0], i[1], i[2])
                    i.append(p)

                #Flatten list
                output_flattened = [item for sublist in nested for item in sublist]

                #Add Accession and GC content information to the list, which will be added to the output CSV
                csv_nested = [[accession], output_flattened, [GC], [GC3]]
                csv_list = [item for sublist in csv_nested for item in sublist]
                csv_total.append(csv_list)

            else:
                continue

    #Set headers for the CSV
    headers = ["Accession", "P0_observed", "Total_UTRs_0", "P0_Expected_freq", "P0_Observed_freq", "P0_Standard_Deviation", "P0_Zscore", "P0__BinomPval",
        "P1_observed", "P1_Total_UTRs", "P1_Expected_freq", "P1_Observed_freq", "P1_Standard_Deviation", "P1_Zscore", "P1__BinomPval",
        "P2_observed", "P2_Total_UTRs", "P2_Expected_freq", "P2_Observed_freq", "P2_Standard_Deviation", "P2_Zscore", "P2__BinomPval",
        "P3_observed", "P3_Total_UTRs", "P3_Expected_freq", "P3_Observed_freq", "P3_Standard_Deviation", "P3_Zscore", "P3__BinomPval",
        "P4_observed", "P4_Total_UTRs", "P4_Expected_freq", "P4_Observed_freq", "P4_Standard_Deviation", "P4_Zscore", "P4__BinomPval",
        "P5_observed", "P5_Total_UTRs", "P5_Expected_freq", "P5_Observed_freq", "P5_Standard_Deviation", "P5_Zscore", "P5__BinomPval",
        "P6_observed", "P6_Total_UTRs", "P6_Expected_freq", "P6_Observed_freq", "P6_Standard_Deviation", "P6_Zscore", "P6__BinomPval",
        "Genomic_GC", "Genomic_GC3"]

    create_csv(headers, csv_total)


def get_utr_stuff(list):

    ''' Obtain the UTR sequence for each gene from the FASTA file '''

    utr_list = []
    total_utr = ''

    #For each gene in the FASTA file...
    for i in list:

        #First generate UTR list, subnested for each codon
        a = i.split("\n")
        utr_seq = a[1]
        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        utr_list.append(codon_seq)

        #Now create total UTR list which will be useful later
        total_utr += utr_seq

    return utr_list, total_utr


def get_freq(utr_list):

    ''' Calculate ASC frequencies at each position to +6 '''

    codons = [0, 0, 0, 0, 0, 0, 0]
    total_utrs = len(utr_list)

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
    frequencies_list = frequencies_array.tolist()

    return codons, total_utrs, frequencies_array, frequencies_list


def get_simulations(utr_list, total_utr, start_time):

    ''' Simulate 10,000 UTR sequences for each genome based upon dinucleotide content of total UTR regions '''

    #Probability of first base - dictionary
    nt_dict = {}
    length = len(total_utr)

    for i in range(length):
        base = total_utr[i]
        if base not in nt_dict:
            nt_dict[base] = 1
        else:
            nt_dict[base] += 1

    nt_freq = {k: v / length for k, v in nt_dict.items()}

    #Second base / Next base
    trans = {}
    for i in range(len(total_utr)-1):

        dinuc = total_utr[i:i+2]
        first_dinuc = dinuc[0]
        sec_dinuc = dinuc[1]

        if first_dinuc not in trans:
            trans[first_dinuc] = [sec_dinuc]
        else:
            trans[first_dinuc] += sec_dinuc

    #Generate simulations
    total_list = []
    count = 0

    #This while loop determines how many simulations will be completed
    while count < 10000:

        sim_n = []
        codon_list = []

        sim = []

        np.random.seed()

        #Calculate next base
        first_base = np.random.choice(['a', 'c', 'g', 't'], p=[nt_freq['a'], nt_freq['c'], nt_freq['g'], nt_freq['t']])
        sim.append(first_base)

        #For one gene simulation
        while (len(sim) < 21):

            #Generate next base
            prev_base = sim[-1]
            next_base = np.random.choice(trans[prev_base], replace=True)
            sim.append(next_base)

        sim = "".join(sim)
        codon_seq = [sim[i:i+3] for i in range(0, len(sim), 3)]
        total_list.append(codon_seq)
        count += 1

    x, y, n_array, n_list = get_freq(total_list)
    end_time = time.time()
    total_time = end_time - start_time
    print ('Simulation complete: ' + str(total_time))

    return n_array, n_list


def create_csv(headers, csv_total):

    ''' Output to CSV file '''

    filename = "output_binom_11_v2_repeat.csv"
    subdir = "4_Outputs"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)


def main():
    folder_parse(source)


### RUN ###

if __name__ == '__main__':
    main()
