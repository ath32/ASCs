### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_mollicutes'

### FUNCTIONS ###

def main(source):

    csv_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            path = os.path.join(source, f)
            raw = open(path).read()

            if '>' in raw:

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
                    a2 = a1[2].strip()
                    a3 = a2.split("=")
                    GC = float(a3[1])

                    a4 = a1[3].strip()
                    a5 = a4.split("=")
                    GC3 = float(a5[1])

                    #Split into sections for use in functions
                    genes = raw.strip().split('>')
                    genes_c = list(filter(None, genes))

                    #Check translation transl_table
                    transl_table = a1[1]

                    if transl_table == ' 4':

                        #Split UTRs into codons
                        glist = get_list(genes_c)

                        #Calculate frequencies for each group of genes
                        ttt_freq_a, ttt_freq_l = get_freqs(glist, 'ttt')
                        ttc_freq_a, ttc_freq_l = get_freqs(glist, 'ttc')
                        tta_freq_a, tta_freq_l = get_freqs(glist, 'tta')
                        ttg_freq_a, ttg_freq_l = get_freqs(glist, 'ttg')

                        tct_freq_a, tct_freq_l = get_freqs(glist, 'tct')
                        tcc_freq_a, tcc_freq_l = get_freqs(glist, 'tcc')
                        tca_freq_a, tca_freq_l = get_freqs(glist, 'tca')
                        tcg_freq_a, tcg_freq_l = get_freqs(glist, 'tcg')

                        tat_freq_a, tat_freq_l = get_freqs(glist, 'tat')
                        tac_freq_a, tac_freq_l = get_freqs(glist, 'tac')
                        taa_freq_a, taa_freq_l = get_freqs(glist, 'taa')
                        tag_freq_a, tag_freq_l = get_freqs(glist, 'tag')

                        tgt_freq_a, tgt_freq_l = get_freqs(glist, 'tgt')
                        tgc_freq_a, tgc_freq_l = get_freqs(glist, 'tgc')
                        tga_freq_a, tga_freq_l = get_freqs(glist, 'tga')
                        tgg_freq_a, tgg_freq_l = get_freqs(glist, 'tgg')

                        ctt_freq_a, ctt_freq_l = get_freqs(glist, 'ctt')
                        ctc_freq_a, ctc_freq_l = get_freqs(glist, 'ctc')
                        cta_freq_a, cta_freq_l = get_freqs(glist, 'cta')
                        ctg_freq_a, ctg_freq_l = get_freqs(glist, 'ctg')

                        cct_freq_a, cct_freq_l = get_freqs(glist, 'cct')
                        ccc_freq_a, ccc_freq_l = get_freqs(glist, 'ccc')
                        cca_freq_a, cca_freq_l = get_freqs(glist, 'cca')
                        ccg_freq_a, ccg_freq_l = get_freqs(glist, 'ccg')

                        cat_freq_a, cat_freq_l = get_freqs(glist, 'cat')
                        cac_freq_a, cac_freq_l = get_freqs(glist, 'cac')
                        caa_freq_a, caa_freq_l = get_freqs(glist, 'caa')
                        cag_freq_a, cag_freq_l = get_freqs(glist, 'cag')

                        cgt_freq_a, cgt_freq_l = get_freqs(glist, 'cgt')
                        cgc_freq_a, cgc_freq_l = get_freqs(glist, 'cgc')
                        cga_freq_a, cga_freq_l = get_freqs(glist, 'cga')
                        cgg_freq_a, cgg_freq_l = get_freqs(glist, 'cgg')

                        att_freq_a, att_freq_l = get_freqs(glist, 'att')
                        atc_freq_a, atc_freq_l = get_freqs(glist, 'atc')
                        ata_freq_a, ata_freq_l = get_freqs(glist, 'ata')
                        atg_freq_a, atg_freq_l = get_freqs(glist, 'atg')

                        act_freq_a, act_freq_l = get_freqs(glist, 'act')
                        acc_freq_a, acc_freq_l = get_freqs(glist, 'acc')
                        aca_freq_a, aca_freq_l = get_freqs(glist, 'aca')
                        acg_freq_a, acg_freq_l = get_freqs(glist, 'acg')

                        aat_freq_a, aat_freq_l = get_freqs(glist, 'aat')
                        aac_freq_a, aac_freq_l = get_freqs(glist, 'aac')
                        aaa_freq_a, aaa_freq_l = get_freqs(glist, 'aaa')
                        aag_freq_a, aag_freq_l = get_freqs(glist, 'aag')

                        agt_freq_a, agt_freq_l = get_freqs(glist, 'agt')
                        agc_freq_a, agc_freq_l = get_freqs(glist, 'agc')
                        aga_freq_a, aga_freq_l = get_freqs(glist, 'aga')
                        agg_freq_a, agg_freq_l = get_freqs(glist, 'agg')

                        gtt_freq_a, gtt_freq_l = get_freqs(glist, 'gtt')
                        gtc_freq_a, gtc_freq_l = get_freqs(glist, 'gtc')
                        gta_freq_a, gta_freq_l = get_freqs(glist, 'gta')
                        gtg_freq_a, gtg_freq_l = get_freqs(glist, 'gtg')

                        gct_freq_a, gct_freq_l = get_freqs(glist, 'gct')
                        gcc_freq_a, gcc_freq_l = get_freqs(glist, 'gcc')
                        gca_freq_a, gca_freq_l = get_freqs(glist, 'gca')
                        gcg_freq_a, gcg_freq_l = get_freqs(glist, 'gcg')

                        gat_freq_a, gat_freq_l = get_freqs(glist, 'gat')
                        gac_freq_a, gac_freq_l = get_freqs(glist, 'gac')
                        gaa_freq_a, gaa_freq_l = get_freqs(glist, 'gaa')
                        gag_freq_a, gag_freq_l = get_freqs(glist, 'gag')

                        ggt_freq_a, ggt_freq_l = get_freqs(glist, 'ggt')
                        ggc_freq_a, ggc_freq_l = get_freqs(glist, 'ggc')
                        gga_freq_a, gga_freq_l = get_freqs(glist, 'gga')
                        ggg_freq_a, ggg_freq_l = get_freqs(glist, 'ggg')

                        #Manipulate data for suitability in CSV file
                        raw_list = [ttt_freq_l, ttc_freq_l, tta_freq_l, ttg_freq_l, tct_freq_l, tcc_freq_l, tca_freq_l, tcg_freq_l, tat_freq_l, tac_freq_l, taa_freq_l, tag_freq_l, tgt_freq_l, tgc_freq_l, tga_freq_l, tgg_freq_l,
                                ctt_freq_l, ctc_freq_l, cta_freq_l, ctg_freq_l, cct_freq_l, ccc_freq_l, cca_freq_l, ccg_freq_l, cat_freq_l, cac_freq_l, caa_freq_l, cag_freq_l, cgt_freq_l, cgc_freq_l, cga_freq_l, cgg_freq_l,
                                att_freq_l, atc_freq_l, ata_freq_l, atg_freq_l, act_freq_l, acc_freq_l, aca_freq_l, acg_freq_l, aat_freq_l, aac_freq_l, aaa_freq_l, aag_freq_l, agt_freq_l, agc_freq_l, aga_freq_l, agg_freq_l,
                                gtt_freq_l, gtc_freq_l, gta_freq_l, gtg_freq_l, gct_freq_l, gcc_freq_l, gca_freq_l, gcg_freq_l, gat_freq_l, gac_freq_l, gaa_freq_l, gag_freq_l, ggt_freq_l, ggc_freq_l, gga_freq_l, ggg_freq_l
                        ]

                        #Flatten list
                        output_flattened = [item for sublist in raw_list for item in sublist]

                        #Create master list, which will represents a row of the output CSV
                        csv_nested = [[accession], output_flattened, [GC], [GC3]]
                        csv_list = [item for sublist in csv_nested for item in sublist]
                        csv_total.append(csv_list)

                    else:
                        continue

    #Set headers for the CSV
    codons = ['ttt', 'ttc', 'tta', 'ttg', 'tct', 'tcc', 'tca', 'tcg', 'tat', 'tac', 'taa', 'tag', 'tgt', 'tgc', 'tga', 'tgg',
                'ctt', 'ctc', 'cta', 'ctg', 'cct', 'ccc', 'cca', 'ccg', 'cat', 'cac', 'caa', 'cag', 'cgt', 'cgc', 'cga', 'cgg',
                'att', 'atc', 'ata', 'atg', 'act', 'acc', 'aca', 'acg', 'aat', 'aac', 'aaa', 'aag', 'agt', 'agc', 'aga', 'agg',
                'gtt', 'gtc', 'gta', 'gtg', 'gct', 'gcc', 'gca', 'gcg', 'gat', 'gac', 'gaa', 'gag', 'ggt', 'ggc', 'gga', 'ggg'
        ]

    headers = []

    headers.append('Accession')

    for i in codons:
        headers.append(i + '_position_0')
        headers.append(i + '_position_1')
        headers.append(i + '_position_2')
        headers.append(i + '_position_3')
        headers.append(i + '_position_4')
        headers.append(i + '_position_5')
        headers.append(i + '_position_6')

    headers.append('GC')
    headers.append('GC3')

    create_csv(headers, csv_total)


def get_list(genes_c):

    ''' Split UTR sequence into codons '''

    utr_list = []

    for i in genes_c:
        b = i.split('\n')
        utr_seq = b[1]

        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        utr_list.append(codon_seq)

    return utr_list


def get_freqs(list, codon):

    ''' Get frequencies for a specific codon '''

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

    ''' Output results to CSV file '''

    filename = "3.5_Frequencies_mollicutes_2.csv"
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
