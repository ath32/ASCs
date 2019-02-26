### IMPORTS ###

import os
import re
import numpy as np

### CHOOSE TARGET GENOME SUCH THAT ONLY THE APPROPRIATE FASTA AND PAXDB IS UNHASHED ###

#A. ferrooxidans - CP001219
fasta = './2_FASTA_CDS_paxdb/CP001219.txt'
PaxDb = './3_PaxDb_data/1. A.ferrooxidans - 243159-WHOLE_ORGANISM-integrated.txt'

#B. anthracis - AE017225
# fasta = './2_FASTA_CDS_paxdb/AE017225.txt'
# PaxDb = './3_PaxDb_data/2. B.anthracis - 260799-GPM_201408.txt'

#B. henselae - BX897699
# fasta = './2_FASTA_CDS_paxdb/BX897699.txt'
# PaxDb = './3_PaxDb_data/3. B.henselae - 283166-Bhenselae_Albrethsen_2013.txt'

#B. thetaiotaomicron - AE015928
# fasta = './2_FASTA_CDS_paxdb/AE015928.txt'
# PaxDb = './3_PaxDb_data/4. B.thetaiotaomicron VPI-5482 - 226186-GPM_201408.txt'

#C. jejuni - AL111168
# fasta = './2_FASTA_CDS_paxdb/AL111168.txt'
# PaxDb = './3_PaxDb_data/5. C.jejuni - 192222-Campylobacter_jejuni_SC_biomart_17684_E__1reps.txt'

#D. deserti - CP001114
# fasta = './2_FASTA_CDS_paxdb/CP001114.txt'
# PaxDb = './3_PaxDb_data/6. D.deserti - 546414-DeinococcusDeserti_PRIDE.txt'

#D. vulgaris - AE017285
# fasta = './2_FASTA_CDS_paxdb/AE017285.txt'
# PaxDb = './3_PaxDb_data/7. D.vulgaris - 882-WHOLE_ORGANISM-integrated.txt'

#E. coli - U00096
# fasta = './2_FASTA_CDS_paxdb/U00096.txt'
# PaxDb = './3_PaxDb_data/8. E.coli - 511145-Ecoli_iBAQ_arike_2012.txt'

#H. pylori - AE000511
# fasta = './2_FASTA_CDS_paxdb/AE000511.txt'
# PaxDb = './3_PaxDb_data/9*. H.pylori 26695 - 85962-WHOLE_ORGANISM-integrated.txt'

#L. interrogans - AE016823
# fasta = './2_FASTA_CDS_paxdb/AE016823.txt'
# PaxDb = './3_PaxDb_data/10. L.interrogans - 267671-PA_201308.txt'

#L. lactis - AE005176
# fasta = './2_FASTA_CDS_paxdb/AE005176.txt'
# PaxDb = './3_PaxDb_data/11*. L.lactis - 272623-Lactococcus_lactis_Lahtvee_2014.txt'

#L. pneumophila - AE017354
# fasta = './2_FASTA_CDS_paxdb/AE017354.txt'
# PaxDb = './3_PaxDb_data/12*. L.pneumophila - 272624-GPM_201408.txt'

#M. aeruginosa - AP009552
# fasta = './2_FASTA_CDS_paxdb/AP009552.txt'
# PaxDb = './3_PaxDb_data/13. M.aeruginosa - 449447-WHOLE_ORGANISM-integrated.txt'

#M. pneumoniae - CP002077*
#fasta = './2_FASTA_CDS_paxdb/CP002077.txt'
#PaxDb = './3_PaxDb_data/14. M.pneumoniae - 722438-Mycoplasma_pneumoniae_M129_Kuhner_et_al_Science2009.txt'

#M. tuberculosis - AL123456
# fasta = './2_FASTA_CDS_paxdb/AL123456.txt'
# PaxDb = './3_PaxDb_data/15*. M.tuberculosis - 83332-WHOLE_ORGANISM-integrated.txt'

#N. meningitidis - AE002098
# fasta = './2_FASTA_CDS_paxdb/AE002098.txt'
# PaxDb = './3_PaxDb_data/16*. N.meningitidis - 122586-GPM_201408.txt'

#S. aureus - BA000017
# fasta = './2_FASTA_CDS_paxdb/BA000017.txt'
# PaxDb = './3_PaxDb_data/17*. S.aureus Mu50 - 158878-GPM_201408.txt'

#S. oneidensis - AE014299
# fasta = './2_FASTA_CDS_paxdb/AE014299.txt'
# PaxDb = './3_PaxDb_data/18. S.oneidensis MR-1 - 211586-GPM_201408.txt'

#S. pyogenes - AE004092
# fasta = './2_FASTA_CDS_paxdb/AE004092.txt'
# PaxDb = './3_PaxDb_data/19*. S.pyogenes M1GAS - 160490-StreptococcusPyogenes_M1GAS_PRIDE.txt'

#S. enterica - AE006468
# fasta = './2_FASTA_CDS_paxdb/AE006468.txt'
# PaxDb = './3_PaxDb_data/20*. S.typhimurium - 99287-WHOLE_ORGANISM-integrated.txt'

#Synechocystis sp. - AP012205
# fasta = './2_FASTA_CDS_paxdb/AP012205.txt'
# PaxDb = './3_PaxDb_data/21*. Synechocystis.sp. 6803 - 1148-SynechocystisSp_strain_PRIDE.txt'

#Y. pestis - AL590842
# fasta = './2_FASTA_CDS_paxdb/AL590842.txt'
# PaxDb = './3_PaxDb_data/22*. Y.pestis - 214092-Yersinia_pestis_SC_biomart_18099_O__3reps.txt'

### FUNCTIONS ###

def main(fasta, PaxDb):

    #Get file paths and load data
    raw_fasta = open(os.path.abspath(fasta)).read()
    raw_paxdb = open(os.path.abspath(PaxDb)).read()

    #Obtain Accession
    g = raw_fasta.split('>')
    g1 = g[1].split(';')
    accession = g1[0]

    #Get gene IDs for two groups
    HEGs, LEGs = get_expression(raw_paxdb)

    #If FASTA gene is in either list, append to two new output strings
    HEGs_fasta, LEGs_fasta = get_new_txt(raw_fasta, HEGs, LEGs)

    #Create two new CSVs
    get_CSVs(HEGs_fasta, LEGs_fasta, accession)


def get_expression(raw_paxdb):

    ''' Obtain the IDs of HEGs and LEGs '''

    #Create nested list for pax data, incl locus tag and protein expression
    raw = []

    #Create list, through which we can access the information for each gene
    if 'raw_spectral_count' in raw_paxdb:
        i = raw_paxdb.split('raw_spectral_count\n')
        i1 = i[1].split('\n')
        del i1[-1]

        for i in i1:
            chunked = i.split('\t')
            internal_id = chunked[0]
            external_id = chunked[1]
            abundance = chunked[2]
            raw_count = chunked[3]
            raw.append([internal_id, external_id, float(abundance)])

    else:
        i = raw_paxdb.split('abundance\n')
        i1 = i[1].split('\n')
        del i1[-1]

        for i in i1:
            chunked = i.split('\t')
            internal_id = chunked[0]
            external_id = chunked[1]
            abundance = chunked[2]
            raw.append([internal_id, external_id, float(abundance)])

    srted = sorted(raw, key = lambda x: float(x[2]))

    #Calculate top and bottom quartiles
    abundances = []

    for i in srted:
        abundances.append(i[2])

    array = np.array(abundances)
    top = np.percentile(array, 75)
    bottom = np.percentile(array, 25)

    #Create two lists - HEGs (top 25%) and LEGs (bottom 25%)
    HEGs_values = []
    HEGs_IDs = []

    for i in srted:
        if i[2] >= top:
            HEGs_values.append(i[2])

    for i in srted:
        if i[2] in HEGs_values:
            HEGs_IDs.append(i[1])

    LEGs_values = []
    LEGs_IDs = []

    for i in srted:
        if i[2] <= bottom:
            LEGs_values.append(i[2])

    for i in srted:
        if i[2] in LEGs_values:
            LEGs_IDs.append(i[1])

    return HEGs_IDs, LEGs_IDs

def get_new_txt(raw_fasta, HEGs, LEGs):

    ''' Build new FASTA files - one for HEGs and one for LEGs for each genome '''

    HEGs_fasta = ''
    LEGs_fasta = ''

    #Refine HEGs and LEGs to just numbers and letters
    refined_HEGs = [''.join(re.findall('[a-zA-Z\d]', i)) for i in HEGs]
    refined_LEGs = [''.join(re.findall('[a-zA-Z\d]', i)) for i in LEGs]

    #Get chunks to iterate through
    chunks = raw_fasta.split('>')
    chunks.pop(0)

    htest = []
    ltest = []

    for i in chunks:
        a = i.split('\n')
        a1 = a[0].split(';')
        id = a1[4]
        refined_id = ''.join(re.findall('[a-zA-Z\d]', id))

        for j in refined_HEGs:
            if refined_id in j:
                HEGs_fasta += ">" + i

        for k in refined_LEGs:
            if refined_id in k:
                LEGs_fasta += ">" + i

    return HEGs_fasta, LEGs_fasta

def get_CSVs(HEGs_fasta, LEGs_fasta, accession):

    ''' Create output files '''

    subdir1 = '2_FASTA_CDS_hegs'
    filename1 = accession + '.txt'
    filepath1 = os.path.join(subdir1, filename1)
    f1 = open(filepath1, 'a')
    f1.write(HEGs_fasta)
    f1.close()

    subdir2 = '2_FASTA_CDS_legs'
    filename2 = accession + '.txt'
    filepath2 = os.path.join(subdir2, filename2)
    f2 = open(filepath2, 'a')
    f2.write(LEGs_fasta)
    f2.close()

### RUN ###

if __name__ == '__main__':
    main(fasta, PaxDb)
