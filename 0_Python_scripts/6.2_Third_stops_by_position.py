### IMPORTS ###

import os
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'
#source = '2_HEGS_fasta'
#source = '2_LEGS_fasta'

### FUNCTIONS ###

def main(source):

    csv_total = []

    for root, dirs, filenames in os.walk(source):

        for f in filenames:
            path = os.path.join(source, f)
            raw = open(path).read()

            #Obtain genome information
            a = raw.split('>')
            a1 = a[1].split(';')

            accession = a1[0]
            gc = float(a1[1].replace('GC', '').replace('=', ''))
            gc3 = float(a1[2].replace('GC', '').replace('=', ''))

            #Access sequences and identify whether there is a stop before position 'n'
            chunks = list(filter(None, a))

            #Split genes into lists for each codon position (can't go beyond 5 because sequences aren't long enough)
            pos1_add, pos1_no_add = get_lists(chunks, 1)
            pos2_add, pos2_no_add = get_lists(chunks, 2)
            pos3_add, pos3_no_add = get_lists(chunks, 3)
            pos4_add, pos4_no_add = get_lists(chunks, 4)
            pos5_add, pos5_no_add = get_lists(chunks, 5)

            #Only use genomes where they possess both types of genes
            if len(pos1_add) < 1 or len(pos2_add) < 1 or len(pos3_add) < 1 or len(pos4_add) < 1 or len(pos5_add) < 1:
                continue
            else:

                #For each set of lists, calculate probability of having a stop at n+1
                add_freq_1, no_add_freq_1, difference_1 = get_third_stop_frequencies(pos1_add, pos1_no_add, 1)
                add_freq_2, no_add_freq_2, difference_2 = get_third_stop_frequencies(pos2_add, pos2_no_add, 2)
                add_freq_3, no_add_freq_3, difference_3 = get_third_stop_frequencies(pos3_add, pos3_no_add, 3)
                add_freq_4, no_add_freq_4, difference_4 = get_third_stop_frequencies(pos4_add, pos4_no_add, 4)
                add_freq_5, no_add_freq_5, difference_5 = get_third_stop_frequencies(pos5_add, pos5_no_add, 5)

                csv_line = [accession, add_freq_1, no_add_freq_1, difference_1, add_freq_2, no_add_freq_2, difference_2, add_freq_3, no_add_freq_3, difference_3, add_freq_4, no_add_freq_4, difference_4, add_freq_5, no_add_freq_5, difference_5]
                csv_total.append(csv_line)

    headers = ["Accession", "Pos1_Stop", "Pos1_No_Stop", "Pos1_difference",
            "Pos2_Stop", "Pos2_No_Stop", "Pos2_difference",
            "Pos3_Stop", "Pos3_No_Stop", "Pos3_difference",
            "Pos4_Stop", "Pos4_No_Stop", "Pos4_difference",
            "Pos5_Stop", "Pos5_No_Stop", "Pos5_difference"
    ]

    create_csv(csv_total, headers)


def get_lists(chunks, position):

    add_stop = []
    no_add_stop = []

    #Edit n for each test
    first_position = 1

    for i in chunks:
        b = i.split('\n')
        utr = b[1]

        #Split into codons
        codon_seq = [utr[i:i+3] for i in range(0, len(utr), 3)]

        #Filter genes into two groups
        if 'taa' in codon_seq[first_position:position+1] or 'tga' in codon_seq[first_position:position+1] or 'tag' in codon_seq[first_position:position+1]:
            add_stop.append(codon_seq)
        else:
            no_add_stop.append(codon_seq)

    return add_stop, no_add_stop


def get_third_stop_frequencies(additional_list, no_additional_list, position):

    add_counter = 0
    no_add_counter = 0

    for i in additional_list:
        if 'taa' in i[position+1:position+2] or 'tga' in i[position+1:position+2] or 'tag' in i[position+1:position+2]:
            add_counter += 1
        else:
            add_counter = add_counter

    for i in no_additional_list:
        if 'taa' in i[position+1:position+2] or 'tga' in i[position+1:position+2] or 'tag' in i[position+1:position+2]:
            no_add_counter += 1
        else:
            no_add_counter = no_add_counter

    #Calculate number of candidate genes for each list
    add_length = len(additional_list)
    no_add_length = len(no_additional_list)

    #Calculate frequency of next stop
    add_freq = add_counter / add_length
    no_add_freq = no_add_counter / no_add_length

    #Calculate difference between the two
    difference = no_add_freq - add_freq

    return add_freq, no_add_freq, difference


def create_csv(csv_total, headers):

    filename = "6.1_Thirdstops_all.csv"
    #filename = "6.1_Thirdstops_HEGs.csv"
    #filename = "6.1_Thirdstops_LEGs.csv"
    subdir = "4_Outputs/CSVs"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main(source)
