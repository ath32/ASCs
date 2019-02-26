### IMPORTS ###

import os
import re
from string import digits
import time
import numpy as np

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_FASTA_Eubacteria_cds_TT11'

def main(source):
    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Define file path
            path = os.path.join(source, f)

            #Define raw data - useful for getting all info from the file
            raw = open(path).read()

            #Create three strings, one for each new txt FASTA file
            taa = ''
            tga = ''
            tag = ''

            #Check file isn't empty
            if '>' in raw:

                #Obtain Accession
                j = raw.split('>')
                j1 = j[1].split(';')
                accession = j1[0].replace(' ', '')

                #Split into sections to access the UTR sequence
                genes = raw.strip().split('>')
                genes_c = list(filter(None, genes))

                #Obtain utrs and based upon their primary stop, build the three fasta files
                for i in genes_c:
                    a = i.split('\n')
                    utr = a[1]
                    primary_stop = utr[:3]

                    if primary_stop == 'taa':
                        taa += ">" + i
                    elif primary_stop == 'tga':
                        tga += ">" + i
                    elif primary_stop == 'tag':
                        tag += ">" + i
                    else:
                        print ('ERROR: no primary stop codon')

            FASTA_creation(taa, tga, tag, accession)


def FASTA_creation(taa, tga, tag, accession):

    ''' Output three unique files based upon the primary stop codon '''

    #Choose a subdirectory, depending on source file used
    subdir1 = '2_FASTA_Eubacteria_cds_TT11_TAA'
    subdir2 = '2_FASTA_Eubacteria_cds_TT11_TGA'
    subdir3 = '2_FASTA_Eubacteria_cds_TT11_TAG'

    filename = accession + '.txt'

    filepath1 = os.path.join(subdir1, filename)
    filepath2 = os.path.join(subdir2, filename)
    filepath3 = os.path.join(subdir3, filename)

    f1 = open(filepath1, 'a')
    f1.write(taa)
    f1.close()

    f2 = open(filepath2, 'a')
    f2.write(tga)
    f2.close()

    f3 = open(filepath3, 'a')
    f3.write(tag)
    f3.close()


### RUN ###

if __name__ == '__main__':
    main(source)
