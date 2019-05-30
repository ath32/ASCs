### IMPORTS ###

import os
import re
from string import digits
import time
import numpy as np
from collections import defaultdict

### CHOOSE SOURCE FOLDER - Unhash the older of interest ###

#ALL
source = 'Aligned'

def main():

    stops = ['taa', 'tga', 'tag']

    n1 = 0
    s1 = 0
    n2 = 0
    s2 = 0
    n3 = 0
    s3 = 0
    n4 = 0
    s4 = 0
    n5 = 0
    s5 = 0
    n6 = 0
    s6 = 0

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Get aligned file (fasta)
            path = os.path.join(source, f)
            raw_file = open(path).read()

            seq_list = []

            split = raw_file.split('>')
            for i in split:
                sections = i.split('\n')
                seq_sections = sections[1:]
                seq = ''.join(seq_sections)
                seq_list.append(seq)

            if seq_list[1] != '' and seq_list[2] != '' and seq_list[3] != '':

                #Get UTRs for all three species (A and B = ingroup species, C = outgroup)
                a, b, c = get_utrs(raw_file)

                #Discard any genes with gaps
                if '-' not in a and '-' not in b and '-' not in c:
                    a_codons = [a[i:i+3] for i in range(0, len(a), 3)]
                    b_codons = [b[i:i+3] for i in range(0, len(b), 3)]
                    c_codons = [c[i:i+3] for i in range(0, len(c), 3)]

                    print (a_codons)

                    #pos1
                    a1 = a_codons[0]
                    b1 = b_codons[0]
                    c1 = c_codons[0]

                    #pos2
                    a2 = a_codons[1]
                    b2 = b_codons[1]
                    c2 = c_codons[1]

                    #pos3
                    a3 = a_codons[2]
                    b3 = b_codons[2]
                    c3 = c_codons[2]

                    #pos4
                    a4 = a_codons[3]
                    b4 = b_codons[3]
                    c4 = c_codons[3]

                    #pos5
                    a5 = a_codons[4]
                    b5 = b_codons[4]
                    c5 = c_codons[4]

                    #pos6
                    a6 = a_codons[5]
                    b6 = b_codons[5]
                    c6 = c_codons[5]

                    #Get switches

                    #Consider all codons of C which are a non-stop - position 1
                    if c1 in stops:
                        #Discard codons where C does not match at least A or B
                        if c1 != a1 and c1 != b1:
                            continue
                            #If C = B = A, we have an ancestral sequence but no codon switch
                        elif c1 == a1 and c1 == b1:
                            n1 += 1
                            s1 += 0
                        #If C = A/B != A/b, we have an ancestral sequence and a switch, which only counts if A/B is not a stop
                        elif c1 == a1 and c1 != b1 and b1 not in stops:
                            n1 += 1
                            s1 += 1
                        elif c1 != a1 and c1 == b1 and a1 not in stops:
                            n1 += 1
                            s1 += 1

                    #Consider all codons of C which are a non-stop - position 2
                    if c2 in stops:
                        #Discard codons where C does not match at least A or B
                        if c2 != a2 and c2 != b1:
                            continue
                            #If C = B = A, we have an ancestral sequence but no codon switch
                        elif c2 == a2 and c2 == b2:
                            n2 += 1
                            s2 += 0
                        #If C = A/B != A/b, we have an ancestral sequence and a switch, which only counts if A/B is a stop
                        elif c2 == a2 and c2 != b2 and b2 not in stops:
                            n2 += 1
                            s2 += 1
                        elif c2 != a2 and c2 == b2 and a2 not in stops:
                            n2 += 1
                            s2 += 1

                    #Consider all codons of C which are a non-stop - position 3
                    if c3 in stops:
                        #Discard codons where C does not match at least A or B
                        if c3 != a3 and c3 != b3:
                            continue
                            #If C = B = A, we have an ancestral sequence but no codon switch
                        elif c3 == a3 and c3 == b3:
                            n3 += 1
                            s3 += 0
                        #If C = A/B != A/b, we have an ancestral sequence and a switch, which only counts if A/B is a stop
                        elif c3 == a3 and c3 != b3 and b3 not in stops:
                            n3 += 1
                            s3 += 1
                        elif c3 != a3 and c3 == b3 and a3 not in stops:
                            n3 += 1
                            s3 += 1

                    #Consider all codons of C which are a non-stop - position 4
                    if c4 in stops:
                        #Discard codons where C does not match at least A or B
                        if c4 != a4 and c4 != b4:
                            continue
                            #If C = B = A, we have an ancestral sequence but no codon switch
                        elif c4 == a4 and c4 == b4:
                            n4 += 1
                            s4 += 0
                        #If C = A/B != A/b, we have an ancestral sequence and a switch, which only counts if A/B is a stop
                        elif c4 == a4 and c4 != b4 and b4 not in stops:
                            n4 += 1
                            s4 += 1
                        elif c4 != a4 and c4 == b4 and a4 not in stops:
                            n4 += 1
                            s4 += 1

                    #Consider all codons of C which are a non-stop - position 5
                    if c5 in stops:
                        #Discard codons where C does not match at least A or B
                        if c5 != a5 and c5 != b5:
                            continue
                            #If C = B = A, we have an ancestral sequence but no codon switch
                        elif c5 == a5 and c5 == b5:
                            n5 += 1
                            s5 += 0
                        #If C = A/B != A/b, we have an ancestral sequence and a switch, which only counts if A/B is a stop
                        elif c5 == a5 and c5 != b5 and b5 not in stops:
                            n5 += 1
                            s5 += 1
                        elif c5 != a5 and c5 == b5 and a5 not in stops:
                            n5 += 1
                            s5 += 1

                    #Consider all codons of C which are a non-stop - position 6
                    if c6 in stops:
                        #Discard codons where C does not match at least A or B
                        if c6 != a6 and c6 != b6:
                            continue
                            #If C = B = A, we have an ancestral sequence but no codon switch
                        elif c6 == a6 and c6 == b6:
                            n6 += 1
                            s6 += 0
                        #If C = A/B != A/b, we have an ancestral sequence and a switch, which only counts if A/B is a stop
                        elif c6 == a6 and c6 != b6 and b6 not in stops:
                            n6 += 1
                            s6 += 1
                        elif c6 != a6 and c6 == b6 and a6 not in stops:
                            n6 += 1
                            s6 += 1

    print (s1, n1, s1/n1)
    print (s2, n2, s2/n2)
    print (s3, n3, s3/n3)
    print (s4, n4, s4/n4)
    print (s5, n5, s5/n5)
    print (s6, n6, s6/n6)

def get_utrs(raw):

    split = raw.split('>')
    utrs = []

    for i in split:
        if i != '':
            chunk = i.split('\n')
            split_seqs = chunk[1:]
            sequence = ''.join(split_seqs)
            if sequence != '':
                utr = sequence[len(sequence)-18:]
                utrs.append(utr)

    if len(utrs) == 3:
        a = utrs[0]
        b = utrs[1]
        c = utrs[2]

    return a, b, c


### RUN ###

if __name__ == '__main__':
    main()
