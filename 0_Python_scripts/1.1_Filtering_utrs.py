### IMPORTS ###

import os
import re
from string import digits
import time
import numpy as np

### CHOOSE SOURCE FOLDER (EMBL) - Unhash the EMBL folder of interest ###

source = '1_Eubacteria_filtered_genomes'
#source = '1_Eubacteria_TT4_genomes'
#source = '3_PaxDb_genomes'

### FUNCTIONS ###

def main(source):
    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Start timer
            start_time = time.time()

            #Define file path
            path = os.path.join(source, f)

            #Define raw data
            raw = open(path).read()

            #Define sequence data - lowercase and no whitespace for accessing sequence data
            s = raw.lower()
            sequence = s.replace(' ', '')

            #Check translation table
            z = raw.split("FT                   /transl_table=")
            z1 = z[1].split("\n")
            transl_table = z1[0]

            #Decide whether to filter for TT11 or TT4 and filter to only use genomes that have 'FT CDS'
            if transl_table == "11":
                if "FT   CDS" in raw:

                    #Obtain Accession
                    a = raw.split('AC')
                    a1 = a[1].split(';')
                    accession = a1[0].replace(' ', '')

                    #Obtain cleaned coding sequence for use in following functions
                    CDS = raw.strip().split('SQ   ')
                    CDS_clean = CDS[1].split('other;')
                    CDS_clean2 = str(CDS_clean[1])
                    CDS_clean3 = CDS_clean2.replace('\n', '').replace(' ', '').replace('//', '')
                    CDS_cleaned = CDS_clean3.translate({ord(k): None for k in digits})

                    #Calculate GC
                    GC = ((CDS_cleaned.count('g') + CDS_cleaned.count('c')) / len(CDS_cleaned) * 100)

                    #Filter genes for those with an IGR > 36 nts
                    gene_samples = IGR_filter(raw, sequence, CDS_cleaned)

                    #Filter genes for CDS checks
                    final_gene_samples, GC3 = CDS_filter(gene_samples, raw, sequence, CDS_cleaned)

                    #Of remaining genes, extract data necessary and export to fasta file
                    get_utr_fasta(final_gene_samples, raw, sequence, CDS_cleaned, accession, GC, GC3)

                    #Obtain time taken per genome
                    stop_time = time.time()
                    total_time = stop_time - start_time
                    print (total_time)


def IGR_filter(raw, sequence, CDS_cleaned):

    ''' Obtain genes which have an IGR of at least 36 nts '''

    #Create lists to be used in the function
    coords_list = []
    gene_samples = []
    final_gene_samples = []

    #Split such that I can access the information from each gene
    i = raw.split('FT   CDS')
    w = len(i)-1

    for t in range(1, w):

        #Obtain locus tag or assign 'locus tag' as a unique identifier
        if 'FT                   /locus_tag="' in i[t]:
            a = i[t].split('FT                   /locus_tag="')
            a1 = a[1].split('"')
            locus_tag = a1[0]
        elif 'FT                   /db_xref="' in i[t]:
            c = i[t].split('FT                   /db_xref="')
            c1 = c[1].split('"')
            locus_tag = c1[0]
        else:
            locus_tag = str(np.random.uniform(-1, 1, size=1))

        #Obtain gene coordinates
        cds_coords = re.findall('\d+', i[t])
        cds_start = int(cds_coords[0])
        cds_stop = int(cds_coords[1])

        #Obtain strand
        strand_string = i[t].replace(" ", "")
        strand = strand_string[0:10]

        if strand == "complement":
            complement = 'complement'
        else:
            complement = 'sense'

        coords_list.append([locus_tag, cds_start, cds_stop, complement])

    #For each gene, calculate the size of the 3' intergenic space
    for i,n in enumerate(coords_list):
        if coords_list[i][3] == 'sense':
            if i != len(coords_list)-1:
                IGR = coords_list[i+1][1] - coords_list[i][2]
                if IGR >= 36:
                    gene_samples.append(coords_list[i][0])

        elif coords_list[i][3] == 'complement':
            if i != 0:
                IGR = coords_list[i][1] - coords_list[i-1][2]
                if IGR >= 36:
                    gene_samples.append(coords_list[i][0])

    return gene_samples


def CDS_filter(gene_samples, raw, sequence, CDS_cleaned):

    ''' Filter genes to ensure the terminate with a stop codon, are length 3n and have no in-frame stops '''

    final_gene_samples = []
    cds_list = []

    #Split such that I can access the information from each gene
    i = raw.split('FT   CDS')
    w = len(i)-1

    for t in range(1, w):

        #Obtain locus tag or assign 'locus tag' as a unique identifier
        if 'FT                   /locus_tag="' in i[t]:
            a = i[t].split('FT                   /locus_tag="')
            a1 = a[1].split('"')
            locus_tag = a1[0]
        elif 'FT                   /db_xref="' in i[t]:
            c = i[t].split('FT                   /db_xref="')
            c1 = c[1].split('"')
            locus_tag = c1[0]
        else:
            locus_tag = str(np.random.uniform(-1, 1, size=1))

        #Obtain gene coordinates
        cds_coords = re.findall('\d+', i[t])
        cds_start = int(cds_coords[0])
        cds_stop = int(cds_coords[1])

        if locus_tag in gene_samples:

            #Obtain CDS
            cds = CDS(sequence, CDS_cleaned, cds_start, cds_stop)

            #Filter out inappropriate gene coding sequences
            if cds_stop_check(cds) == False:
                continue
            elif length_check(cds) == False:
                continue
            elif check_inframe_stop(cds) == False:
                continue
            else:
                final_gene_samples.append(locus_tag)
                cds_list.append(cds)

    #Calclate GC3 content of each gene
    gc3_values = []

    for i in cds_list:

        #Get every third base of each gene cds
        gc3_string = i[2::3]

        #Calculate GC content of this new string
        gc3_cds = ((gc3_string.count('g') + gc3_string.count('c')) / len(gc3_string) * 100)
        gc3_values.append(gc3_cds)

    GC3 = np.mean(gc3_values)

    return final_gene_samples, GC3


def get_utr_fasta(final_gene_samples, raw, sequence, CDS_cleaned, accession, GC, GC3):

    ''' Obtain the 3' UTR sequence (length 30 nt including stop codon) '''

    #Split such that I can access the information from each gene
    i = raw.split('FT   CDS')
    w = len(i)-1

    #Create FASTA string to be added to later
    fasta = ""

    for t in range(1, w):

        #Obtain locus tag, or assign gene id
        if 'FT                   /locus_tag="' in i[t]:
            a = i[t].split('FT                   /locus_tag="')
            a1 = a[1].split('"')
            locus_tag = a1[0]
        elif 'FT                   /db_xref="' in i[t]:
            c = i[t].split('FT                   /db_xref="')
            c1 = c[1].split('"')
            locus_tag = c1[0]
        else:
            locus_tag = str(np.random.uniform(-1, 1, size=1))

        #Only use genes that have come through IGR and CDS filtering
        if locus_tag in final_gene_samples:

            #Obtain gene id
            if 'FT                   /gene="' in i[t]:
                g = i[t].split('FT                   /gene="')
                g1 = g[1].split('"')
                gene_id = g1[0]
            else:
                gene_id = 'Unknown_GeneID'

            #Obtain gene coordinates
            cds_coords = re.findall('\d+', i[t])
            cds_start = int(cds_coords[0])
            cds_stop = int(cds_coords[1])

            #Obtain locus tag
            if 'FT                   /locus_tag="' in i[t]:
                a = i[t].split('FT                   /locus_tag="')
                a1 = a[1].split('"')
                locus_tag = a1[0]
            else:
                locus_tag = 'Unknown_LocusID'

            #Obtain protein id
            if 'FT                   /protein_id="' in i[t]:
                c = i[t].split('FT                   /protein_id="')
                c1 = c[1].split('"')
                protein_id = c1[0]
            else:
                protein_id = "Unknown_ProteinID"

            #Obtain 3' UTR nucleotide sequence
            utr_seq = utr_sequence(raw, sequence, cds_start, cds_stop, CDS_cleaned)

            #Write FASTA_section
            FASTA_section = ">" + accession + "; " + "GC=" + str(GC) + "; GC3=" + str(GC3) + "; " + gene_id + "; " + locus_tag + "; "+ protein_id + ";\n" + utr_seq + "\n"

            #Check UTR sequences are suitable for analysis
            if utr_seq == None:
                continue
            elif start_check(utr_seq) == False:
                continue
            elif nucleotide_check(utr_seq) == False:
                continue
            elif length_check(utr_seq) == False:
                continue
            else:
                fasta += FASTA_section

    #Create output FASTA file
    FASTA_creation(fasta, accession)


def CDS(sequence, CDS_cleaned, cds_start, cds_stop):

    ''' Obtain coding sequence, given start and stop coordinates '''

    n = str(cds_start)
    m = str(cds_stop)

    i = sequence.strip().split('(' + n + '..' + m + ')')
    l = len(i[0])

    for w in i:
        if w[l-10:l] == 'complement':
            rc = CDS_cleaned[cds_start-1:cds_stop].replace('a', 'T').replace('t', 'A').replace('c', 'G').replace('g', 'C').lower()
            reverse_complement = rc[::-1]
            return reverse_complement
        else:
            cds = CDS_cleaned[cds_start-1:cds_stop]
            return cds


def utr_sequence(raw, sequence, cds_start, cds_stop, CDS_cleaned):

    ''' Obtain UTR sequence, given start and stop coordinates '''

    n = str(cds_start)
    m = str(cds_stop)

    sense_utr_start = cds_stop - 2
    sense_utr_stop = sense_utr_start + 29

    antisense_utr_stop = cds_start + 2
    antisense_utr_start = antisense_utr_stop - 29

    i = sequence.split('(' + n + '..' + m + ')')
    l = len(i[0])

    for w in i:
        if w[l-10:l] == 'complement':
            complement = CDS_cleaned[antisense_utr_start-1:antisense_utr_stop].replace('a', 'T').replace('t', 'A').replace('c', 'G').replace('g', 'C').lower()
            reverse_complement = complement[::-1]
            return reverse_complement

        else:
            utr = CDS_cleaned[sense_utr_start-1:sense_utr_stop]
            return utr


def FASTA_creation(fasta, accession):

    ''' Create FASTA output file '''

    #Choose a subdirectory, depending on source file used
    subdir = '2_FASTA_Eubacteria_cds_TT11'
    #subdir = '2_FASTA_Eubacteria_cds_TT4'
    #subdir = '2_FASTA_PaxDb'
    filename = accession + '.txt'
    filepath = os.path.join(subdir, filename)

    f = open(filepath, 'a')

    f.write(fasta)
    f.close()

### UTR CHECK FUNCTIONS ###

def length_check(sequence):
    t = len(sequence)
    if t % 3 == 0:
        length = True
    else:
        length = False
    return length

def start_check(sequence):
    if sequence[0:3] == 'tag':
        start = True
    elif sequence[0:3] == 'tga':
        start = True
    elif sequence[0:3] == 'taa':
        start = True
    else:
        start = False
    return start

def nucleotide_check(sequence):
    filtered = re.findall('[actg]', sequence)
    f_string = "".join(filtered)

    if f_string == sequence:
        nucleotides = True
    else:
        nucleotides = False

    return nucleotides

def cds_stop_check(sequence):
    s = len(sequence)
    if sequence[s-3:s] == 'tag':
        stop = True
    elif sequence[s-3:s] == 'taa':
        stop = True
    elif sequence[s-3:s] == 'tga':
        stop = True
    else:
        stop = False
    return stop

def check_inframe_stop(sequence):

	filter_check = True
	stop_codons = {}
	stop_codons[11] = ['taa', 'tga', 'tag']

	for i in range(0, len(sequence)-3, 3):
		if sequence[i:i+3] in stop_codons[11]:
			filter_check = False

	return filter_check


### RUN ###

if __name__ == '__main__':
    main(source)
