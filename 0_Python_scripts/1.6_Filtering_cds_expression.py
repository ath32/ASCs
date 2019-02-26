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
# fasta = './2_FASTA_CDS_paxdb/CP002077.txt'
# PaxDb = './3_PaxDb_data/14. M.pneumoniae - 722438-Mycoplasma_pneumoniae_M129_Kuhner_et_al_Science2009.txt'

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
#fasta = './2_FASTA_CDS_paxdb/AL590842.txt'
#PaxDb = './3_PaxDb_data/22*. Y.pestis - 214092-Yersinia_pestis_SC_biomart_18099_O__3reps.txt'

### FUNCTIONS ###

def main(fasta, PaxDb):

    #Get file paths and load data
    raw_fasta = open(os.path.abspath(fasta)).read()
    raw_paxdb = open(os.path.abspath(PaxDb)).read()

    #Obtain Accession
    g = raw_fasta.split('>')
    g1 = g[1].split(';')
    accession = g1[0]

    #Get ID lists
    paxdb_ids = get_paxdb_id(raw_paxdb)

    #Get FASTA
    FASTA_total = get_fasta(paxdb_ids, raw_fasta)

    #Make new file
    get_FASTA(FASTA_total, accession)


def get_paxdb_id(raw_paxdb):

    list = []

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
            external_id_refined = ''.join(re.findall('[a-zA-Z\d]', external_id))
            list.append([internal_id, external_id_refined, float(abundance)])

    else:
        i = raw_paxdb.split('abundance\n')
        i1 = i[1].split('\n')
        del i1[-1]

        for i in i1:
            chunked = i.split('\t')
            internal_id = chunked[0]
            external_id = chunked[1]
            abundance = chunked[2]
            external_id_refined = ''.join(re.findall('[a-zA-Z\d]', external_id))
            list.append([internal_id, external_id_refined, float(abundance)])

    return list


def get_fasta(paxdb_ids, raw_fasta):

    FASTA_total = ''

    #Get chunks to iterate through
    chunks = raw_fasta.split('>')
    chunks.pop(0)

    for i in chunks:

        #Get FASTA ID and FASTA chunk
        a = i.split('\n')
        a1 = a[0].split(';')
        id = a1[4]
        refined_id = ''.join(re.findall('[a-zA-Z\d]', id))

        #Fasta line
        fasta_line = a[0]
        sequence = a[1]

        for i in paxdb_ids:

            #Get PaxDb info
            external_id = i[1]
            abundance = i[2]

            if refined_id in external_id:

                new_fasta = ">" + fasta_line + " " + str(abundance) + ";\n" + sequence + "\n"
                FASTA_total += new_fasta

    return FASTA_total


def get_FASTA(FASTA_total, accession):

    ''' Create output files '''

    subdir = '2_FASTA_CDS_expression'
    filename = accession + '.txt'
    filepath = os.path.join(subdir, filename)
    f1 = open(filepath, 'a')
    f1.write(FASTA_total)
    f1.close()


### RUN ###

if __name__ == '__main__':
    main(fasta, PaxDb)
