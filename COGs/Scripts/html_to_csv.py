### IMPORTS ###

import os
import re
import csv

#SELECT COG

#ATGC001
# html1 = os.path.abspath("ATGC001/COG_membership/html/Salmonella_bongori_NCTC_12419.GCF_000252995.1.html")
# html2 = os.path.abspath("ATGC001/COG_membership/html/Salmonella_enterica_arizonae_serovar_62_z36_str_RKS2983.GCF_000756465.1.html")
# html3 = os.path.abspath("ATGC001/COG_membership/html/Salmonella_enterica_enterica_serovar_Choleraesuis.GCF_000742815.1.html")

#ATGC003
# html1 = os.path.abspath("ATGC003/COG_membership/html/Streptococcus_pneumoniae_INV200.GCF_000210935.1.html")
# html2 = os.path.abspath("ATGC003/COG_membership/html/Streptococcus_pseudopneumoniae_IS7493.GCF_000221985.1.html")
# html3 = os.path.abspath("ATGC003/COG_membership/html/Streptococcus_mitis_B6.GCF_000027165.1.html")

#ATGC008
# html1 = os.path.abspath("ATGC008/COG_membership/html/Streptococcus_salivarius_CCHSS3.GCF_000253335.1.html")
# html2 = os.path.abspath("ATGC008/COG_membership/html/Streptococcus_salivarius_JIM8777.GCF_000253315.1.html")
# html3 = os.path.abspath("ATGC008/COG_membership/html/Streptococcus_thermophilus_CNRZ1066.GCF_000011845.1.html")

#ATGC014 - ready
# html1 = os.path.abspath("ATGC014/COG_membership/html/Bacillus_cereus_ATCC_10987.GCF_000008005.1.html")
# html2 = os.path.abspath("ATGC014/COG_membership/html/Bacillus_thuringiensis_Al_Hakam.GCF_000015065.1.html")
# html3 = os.path.abspath("ATGC014/COG_membership/html/Bacillus_weihenstephanensis_KBAB4.GCF_000018825.1.html")

#ATGC015 - ready
# html1 = os.path.abspath("ATGC015/COG_membership/html/Bacillus_JS.GCF_000259365.1.html")
# html2 = os.path.abspath("ATGC015/COG_membership/html/Bacillus_subtilis_subtilis_str_BSP1.GCF_000321395.1.html")
# html3 = os.path.abspath("ATGC015/COG_membership/html/Bacillus_subtilis_spizizenii_str_W23.GCF_000146565.1.html")

#ATGC035 - ready
# html1 = os.path.abspath("ATGC035/COG_membership/html/Mycoplasma_capricolum_capricolum_ATCC_27343.GCF_000012765.1.html")
# html2 = os.path.abspath("ATGC035/COG_membership/html/Mycoplasma_leachii_PG50.GCF_000183365.1.html")
# html3 = os.path.abspath("ATGC035/COG_membership/html/Mycoplasma_mycoides_mycoides_SC_str_PG1.GCF_000011445.1.html")

#ATGC044 - ready
# html1 = os.path.abspath("ATGC044/COG_membership/html/Rickettsia_akari_Hartford.GCF_000018205.1.html")
# html2 = os.path.abspath("ATGC044/COG_membership/html/Rickettsia_massiliae_MTU5.GCF_000016625.1.html")
# html3 = os.path.abspath("ATGC044/COG_membership/html/Rickettsia_prowazekii_Rp22.GCF_000022785.1.html")

#ATGC050 - ready
# html1 = os.path.abspath("ATGC050/COG_membership/html/Helicobacter_pylori_83.GCF_000213135.1.html")
# html2 = os.path.abspath("ATGC050/COG_membership/html/Helicobacter_pylori_SouthAfrica7.GCF_000185245.1.html")
# html3 = os.path.abspath("ATGC050/COG_membership/html/Helicobacter_cetorum_MIT_99_5656.GCF_000259275.1.html")

#ATGC071 - ready
# html1 = os.path.abspath("ATGC071/COG_membership/html/Pseudomonas_putida_H8234.GCF_000410575.1.html")
# html2 = os.path.abspath("ATGC071/COG_membership/html/Pseudomonas_putida_NBRC_14164.GCF_000412675.1.html")
# html3 = os.path.abspath("ATGC071/COG_membership/html/Pseudomonas_putida_GB_1.GCF_000019125.1.html")

#ATGC088 - ready
# html1 = os.path.abspath("ATGC088/COG_membership/html/Burkholderia_thailandensis_E264.GCF_000012365.1.html")
# html2 = os.path.abspath("ATGC088/COG_membership/html/Burkholderia_thailandensis_MSMB121.GCF_000385525.1.html")
# html3 = os.path.abspath("ATGC088/COG_membership/html/Burkholderia_pseudomallei_K96243.GCF_000011545.1.html")

#ATGC089 - ready
# html1 = os.path.abspath("ATGC089/COG_membership/html/Burkholderia_cenocepacia_J2315.GCF_000009485.1.html")
# html2 = os.path.abspath("ATGC089/COG_membership/html/Burkholderia_cenocepacia_MC0_3.GCF_000019505.1.html")
# html3 = os.path.abspath("ATGC089/COG_membership/html/Burkholderia_cepacia_GG4.GCF_000292915.1.html")

#ATGC097 - ready
# html1 = os.path.abspath("ATGC097/COG_membership/html/Erwinia_Ejp617.GCF_000165815.1.html")
# html2 = os.path.abspath("ATGC097/COG_membership/html/Erwinia_amylovora_ATCC_49946.GCF_000027205.1.html")
# html3 = os.path.abspath("ATGC097/COG_membership/html/Erwinia_tasmaniensis_Et1_99.GCF_000026185.1.html")

#ATGC100 - ready
# html1 = os.path.abspath("ATGC100/COG_membership/html/Pectobacterium_carotovorum_carotovorum_PC1.GCF_000023605.1.html")
# html2 = os.path.abspath("ATGC100/COG_membership/html/Pectobacterium_carotovorum_carotovorum_PCC21.GCF_000294535.1.html")
# html3 = os.path.abspath("ATGC100/COG_membership/html/Pectobacterium_SCC3193.GCF_000260925.1.html")

#ATGC104 - ready
# html1 = os.path.abspath("ATGC104/COG_membership/html/Bifidobacterium_longum_longum_JCM_1217.GCF_000196555.1.html")
# html2 = os.path.abspath("ATGC104/COG_membership/html/Bifidobacterium_longum_infantis_ATCC_15697_JCM_1222_DSM_20088.GCF_000020425.1.html")
# html3 = os.path.abspath("ATGC104/COG_membership/html/Bifidobacterium_breve_ACS_071_V_Sch8b.GCF_000213865.1.html")

#ATGC108 - can't find outgroup genome in database ******************

#ATGC111 - ready
# html1 = os.path.abspath("ATGC111/COG_membership/html/Aeromonas_hydrophila_hydrophila_ATCC_7966.GCF_000014805.1.html")
# html2 = os.path.abspath("ATGC111/COG_membership/html/Aeromonas_hydrophila_ML09_119.GCF_000401555.1.html")
# html3 = os.path.abspath("ATGC111/COG_membership/html/Aeromonas_salmonicida_salmonicida_A449.GCF_000196395.1.html")

#ATGC123 - insufficient information to differentiate between species 1 and species 3 ***************

#ATGC125 - missing from database ******************

#ATGC134 - ready
# html1 = os.path.abspath("ATGC134/COG_membership/html/Xanthomonas_axonopodis_citri_str_306.GCF_000007165.1.html")
# html2 = os.path.abspath("ATGC134/COG_membership/html/Xanthomonas_axonopodis_citrumelo_F1.GCF_000225915.1.html")
# html3 = os.path.abspath("ATGC134/COG_membership/html/Xanthomonas_oryzae_oryzicola_BLS256.GCF_000168315.3.html")

#ATGC135 - ready
# html1 = os.path.abspath("ATGC135/COG_membership/html/Stenotrophomonas_maltophilia_D457.GCF_000284595.1.html")
# html2 = os.path.abspath("ATGC135/COG_membership/html/Stenotrophomonas_maltophilia_JV3.GCF_000223885.1.html")
# html3 = os.path.abspath("ATGC135/COG_membership/html/Stenotrophomonas_maltophilia_R551_3.GCF_000020665.1.html")

#ATGC137 - species 1 missing from database ********************

#ATGC138 - ready
# html1 = os.path.abspath("ATGC138/COG_membership/html/Francisella_noatunensis_orientalis.GCF_001042565.1.html")
# html2 = os.path.abspath("ATGC138/COG_membership/html/Francisella_philomiragia_philomiragia_ATCC_25017.GCF_000019285.1.html")
# html3 = os.path.abspath("ATGC138/COG_membership/html/Francisella_TX077308.GCF_000219045.1.html")

#ATGC144 - ready
# html1 = os.path.abspath("ATGC144/COG_membership/html/Borreliella_bissettii_DN127.GCF_000222305.1.html")
# html2 = os.path.abspath("ATGC144/COG_membership/html/Borrelia_burgdorferi_ZS7.GCF_000021405.1.html")
# html3 = os.path.abspath("ATGC144/COG_membership/html/Borrelia_garinii_BgVir.GCF_000239475.1.html")

#ATGC147 - ready
# html1 = os.path.abspath("ATGC147/COG_membership/html/Methanococcus_maripaludis_C5.GCF_000016125.1.html")
# html2 = os.path.abspath("ATGC147/COG_membership/html/Methanococcus_maripaludis_C6.GCF_000018485.1.html")
# html3 = os.path.abspath("ATGC147/COG_membership/html/Methanococcus_maripaludis_C7.GCF_000017225.1.html")

#ATGC149 - ready
# html1 = os.path.abspath("ATGC149/COG_membership/html/Acinetobacter_calcoaceticus_PHEA_2.GCF_000191145.1.html")
# html2 = os.path.abspath("ATGC149/COG_membership/html/Acinetobacter_oleivorans_DR1.GCF_000196795.1.html")
# html3 = os.path.abspath("ATGC149/COG_membership/html/Acinetobacter_baumannii_AB307_0294.GCF_000021145.1.html")

#ATGC165 - ready
# html1 = os.path.abspath("ATGC165/COG_membership/html/Rhodobacter_sphaeroides_ATCC_17029.GCF_000015985.1.html")
# html2 = os.path.abspath("ATGC165/COG_membership/html/Rhodobacter_sphaeroides_KD131.GCF_000021005.1.html")
# html3 = os.path.abspath("ATGC165/COG_membership/html/Rhodobacter_sphaeroides_ATCC_17025.GCF_000016405.1.html")

#ATGC171 - ready
# html1 = os.path.abspath("ATGC171/COG_membership/html/Thermoanaerobacter_X514.GCF_000019065.1.html")
# html2 = os.path.abspath("ATGC171/COG_membership/html/Thermoanaerobacter_wiegelii_Rt8_B1.GCF_000147695.2.html")
# html3 = os.path.abspath("ATGC171/COG_membership/html/Thermoanaerobacter_mathranii_mathranii_str_A3.GCF_000092965.1.html")

#ATGC177 - ready
# html1 = os.path.abspath("ATGC177/COG_membership/html/Prochlorococcus_marinus_AS9601.GCF_000015645.1.html")
# html2 = os.path.abspath("ATGC177/COG_membership/html/Prochlorococcus_marinus_MIT_9301.GCF_000015965.1.html")
# html3 = os.path.abspath("ATGC177/COG_membership/html/Prochlorococcus_marinus_MIT_9312.GCF_000012645.1.html")

#ATGC181 - ready
# html1 = os.path.abspath("ATGC181/COG_membership/html/Caldicellulosiruptor_bescii_DSM_6725.GCF_000022325.1.html")
# html2 = os.path.abspath("ATGC181/COG_membership/html/Caldicellulosiruptor_hydrothermalis_108.GCF_000166355.1.html")
# html3 = os.path.abspath("ATGC181/COG_membership/html/Caldicellulosiruptor_owensensis_OL.GCF_000166335.1.html")

#ATGC188 - outgroup is missing from database ***********

#ATGC189 - outgroup is missing from database ***********

#ATGC199 - ready
# html1 = os.path.abspath("ATGC199/COG_membership/html/Paenibacillus_polymyxa_E681.GCF_000146875.3.html")
# html2 = os.path.abspath("ATGC199/COG_membership/html/Paenibacillus_polymyxa_SC2.GCF_000164985.3.html")
# html3 = os.path.abspath("ATGC199/COG_membership/html/Paenibacillus_terrae_HPL_003.GCF_000235585.1.html")

#ATGC201 - ready
# html1 = os.path.abspath("ATGC201/COG_membership/html/Bartonella_grahamii_as4aup.GCF_000022725.1.html")
# html2 = os.path.abspath("ATGC201/COG_membership/html/Bartonella_tribocorum_CIP_105476.GCF_000196435.1.html")
# html3 = os.path.abspath("ATGC201/COG_membership/html/Bartonella_quintana_Toulouse.GCF_000046685.1.html")

#ATGC210 - species 2 is missing from database *************

#ATGC213 - ready
# html1 = os.path.abspath("ATGC213/COG_membership/html/Methylobacterium_extorquens_AM1.GCF_000022685.1.html")
# html2 = os.path.abspath("ATGC213/COG_membership/html/Methylobacterium_extorquens_PA1.GCF_000018845.1.html")
# html3 = os.path.abspath("ATGC213/COG_membership/html/Methylobacterium_populi_BJ001.GCF_000019945.1.html")

#ATGC234 - ready
# html1 = os.path.abspath("ATGC234/COG_membership/html/Clavibacter_michiganensis_michiganensis_NCPPB_382.GCF_000063485.1.html")
# html2 = os.path.abspath("ATGC234/COG_membership/html/Clavibacter_michiganensis_nebraskensis_NCPPB_2581.GCF_000355695.1.html")
# html3 = os.path.abspath("ATGC234/COG_membership/html/Clavibacter_michiganensis_sepedonicus.GCF_000069225.1.html")

#ATGC252 - ready
# html1 = os.path.abspath("ATGC252/COG_membership/html/Anaeromyxobacter_K.GCF_000020805.1.html")
# html2 = os.path.abspath("ATGC252/COG_membership/html/Anaeromyxobacter_dehalogenans_2CP_1.GCF_000022145.1.html")
# html3 = os.path.abspath("ATGC252/COG_membership/html/Anaeromyxobacter_dehalogenans_2CP_C.GCF_000013385.1.html")


def main():
    functions(html1)
    functions(html2)
    functions(html3)

def functions(html_file):

    #Get raw HTML
    raw = open(html_file).read()

    #Split by each row
    split = raw.split('<tr><td class="headcolumn">')

    csv_total = []
    headers = ['Protein_ID', 'COG', 'COG_match_level', 'Synteny', 'Contig_accession', 'Coordinates', 'Protein_length',
    'Locus', 'Protein_accession', 'Product_name', 'Start_codon', 'Stop_codon']

    #Split HTML file into sections required for the analysis
    for i in split:
        if '<html>' not in i:
            j = i.replace('</td><td>', ',').replace('</td></tr>', '').replace('</tbody></table>', '').replace('</body></html>', '')
            list = j.split(',')
            protein_id = list[0]
            cog = list[1]
            cog_match_level = list[2]
            synteny = list[3]
            contig_accession = list[4]
            coords = list[5]
            protein_length = list[6]
            locus = list[7]
            protein_accession = list[8]
            product_name = list[9]
            start = list[13]
            stop = list[14]

            #Retain information necessary for analysis, in CSV format
            row = [protein_id, cog, cog_match_level, synteny, contig_accession, coords, protein_length, locus, protein_accession, product_name, start, stop]
            csv_total.append(row)

    #Get file name and path
    title_split = html_file.split('/')
    name = title_split[9].replace('.html', '')
    folder = title_split[6]

    filename = name + ".csv"
    subdir = folder + "/COG_membership/csv"
    filepath = os.path.join(subdir, filename)

    #Write to CSV file
    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for x in csv_total:
            writer.writerow(x)

### RUN ###

if __name__ == '__main__':
    main()
