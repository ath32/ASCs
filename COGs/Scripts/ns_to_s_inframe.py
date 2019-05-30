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
    n = 0
    s = 0

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Get aligned file (fasta)
            path = os.path.join(source, f)
            raw_file = open(path).read()

            #Get UTRs for all three species (A and B = ingroup species, C = outgroup)
            a, b, c = get_utrs(raw_file)

            #Discard any genes with gaps
            if '-' not in a and '-' not in b and '-' not in c:
                a_codons = [a[i:i+3] for i in range(3, len(a), 3)]
                b_codons = [b[i:i+3] for i in range(3, len(b), 3)]
                c_codons = [c[i:i+3] for i in range(3, len(c), 3)]

                #Consider all codons of C which are a non-stop
                for i,j in enumerate(c_codons):
                    if c_codons[i] not in stops:
                        #Discard codons where C does not match at least A or B
                        if c_codons[i] != a_codons[i] and c_codons[i] != b_codons[i]:
                            continue
                        #If C = B = A, we have an ancestral sequence but no codon switch
                        elif c_codons[i] == a_codons[i] and c_codons[i] == b_codons[i]:
                            n += 1
                            s += 0
                        #If C = A/B != A/b, we have an ancestral sequence and a switch, which only counts if A/B is a stop
                        elif c_codons[i] == a_codons[i] and c_codons[i] != b_codons[i] and b_codons[i] in stops:
                            n += 1
                            s += 1
                        elif c_codons[i] != a_codons[i] and c_codons[i] == b_codons[i] and a_codons[i] in stops:
                            n += 1
                            s += 1
    print (s, n, s/n)

def get_utrs(raw):

    split = raw.split('>')
    utrs = []

    for i in split:
        if i != '':
            chunk = i.split('\n')
            split_seqs = chunk[1:]
            sequence = ''.join(split_seqs)
            utr = sequence[len(sequence)-18:]
            utrs.append(utr)

    a = utrs[0]
    b = utrs[1]
    c = utrs[2]

    return a, b, c


### RUN ###

if __name__ == '__main__':
    main()
