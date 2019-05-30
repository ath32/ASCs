### IMPORTS ###

import os
import re
from string import digits
import time
import numpy as np
from collections import defaultdict

### CHOOSE SOURCE FOLDER (EMBL) - Unhash the EMBL folder of interest ###

#ATGC001
# csv3 = os.path.abspath("ATGC001/COG_membership/csv/Salmonella_bongori_NCTC_12419.GCF_000252995.1.csv")
# csv1 = os.path.abspath("ATGC001/COG_membership/csv/Salmonella_enterica_arizonae_serovar_62_z36_str_RKS2983.GCF_000756465.1.csv")
# csv2 = os.path.abspath("ATGC001/COG_membership/csv/Salmonella_enterica_enterica_serovar_Choleraesuis.GCF_000742815.1.csv")
# wgs3 = os.path.abspath("ATGC001/Genomes/NC_015761.1.fa")
# wgs1 = os.path.abspath("ATGC001/Genomes/NZ_CP006693.1.fa")
# wgs2 = os.path.abspath("ATGC001/Genomes/NZ_CP007639.1.fa")

#ATGC003
# csv1 = os.path.abspath("ATGC003/COG_membership/csv/Streptococcus_pneumoniae_INV200.GCF_000210935.1.csv")
# csv2 = os.path.abspath("ATGC003/COG_membership/csv/Streptococcus_pseudopneumoniae_IS7493.GCF_000221985.1.csv")
# csv3 = os.path.abspath("ATGC003/COG_membership/csv/Streptococcus_mitis_B6.GCF_000027165.1.csv")
# wgs1 = os.path.abspath("ATGC003/Genomes/NC_017593.1.fa")
# wgs2 = os.path.abspath("ATGC003/Genomes/NC_015875.1.fa")
# wgs3 = os.path.abspath("ATGC003/Genomes/NC_013853.1.fa")

#ATGC008
# csv1 = os.path.abspath("ATGC008/COG_membership/csv/Streptococcus_salivarius_CCHSS3.GCF_000253335.1.csv")
# csv2 = os.path.abspath("ATGC008/COG_membership/csv/Streptococcus_salivarius_JIM8777.GCF_000253315.1.csv")
# csv3 = os.path.abspath("ATGC008/COG_membership/csv/Streptococcus_thermophilus_CNRZ1066.GCF_000011845.1.csv")
# wgs1 = os.path.abspath("ATGC008/Genomes/NC_015760.1.fa")
# wgs2 = os.path.abspath("ATGC008/Genomes/NC_017595.1.fa")
# wgs3 = os.path.abspath("ATGC008/Genomes/NC_006449.1.fa")

#ATGC014 - ready
# csv1 = os.path.abspath("ATGC014/COG_membership/csv/Bacillus_cereus_ATCC_10987.GCF_000008005.1.csv")
# csv2 = os.path.abspath("ATGC014/COG_membership/csv/Bacillus_thuringiensis_Al_Hakam.GCF_000015065.1.csv")
# csv3 = os.path.abspath("ATGC014/COG_membership/csv/Bacillus_weihenstephanensis_KBAB4.GCF_000018825.1.csv")
# wgs1 = os.path.abspath("ATGC014/Genomes/NC_003909.8.fa")
# wgs2 = os.path.abspath("ATGC014/Genomes/NC_008600.1.fa")
# wgs3 = os.path.abspath("ATGC014/Genomes/NC_010184.1.fa")

#ATGC015 - ready
# csv1 = os.path.abspath("ATGC015/COG_membership/csv/Bacillus_JS.GCF_000259365.1.csv")
# csv2 = os.path.abspath("ATGC015/COG_membership/csv/Bacillus_subtilis_subtilis_str_BSP1.GCF_000321395.1.csv")
# csv3 = os.path.abspath("ATGC015/COG_membership/csv/Bacillus_subtilis_spizizenii_str_W23.GCF_000146565.1.csv")
# wgs1 = os.path.abspath("ATGC015/Genomes/NC_017743.1.fa")
# wgs2 = os.path.abspath("ATGC015/Genomes/NC_019896.1.fa")
# wgs3 = os.path.abspath("ATGC015/Genomes/NC_014479.1.fa")

#ATGC035 - ready
# csv1 = os.path.abspath("ATGC035/COG_membership/csv/Mycoplasma_capricolum_capricolum_ATCC_27343.GCF_000012765.1.csv")
# csv2 = os.path.abspath("ATGC035/COG_membership/csv/Mycoplasma_leachii_PG50.GCF_000183365.1.csv")
# csv3 = os.path.abspath("ATGC035/COG_membership/csv/Mycoplasma_mycoides_mycoides_SC_str_PG1.GCF_000011445.1.csv")
# wgs1 = os.path.abspath("ATGC035/Genomes/NC_007633.1.fa")
# wgs2 = os.path.abspath("ATGC035/Genomes/NC_014751.1.fa")
# wgs3 = os.path.abspath("ATGC035/Genomes/NC_005364.2.fa")

#ATGC044 - ready
# csv1 = os.path.abspath("ATGC044/COG_membership/csv/Rickettsia_akari_Hartford.GCF_000018205.1.csv")
# csv2 = os.path.abspath("ATGC044/COG_membership/csv/Rickettsia_massiliae_MTU5.GCF_000016625.1.csv")
# csv3 = os.path.abspath("ATGC044/COG_membership/csv/Rickettsia_prowazekii_Rp22.GCF_000022785.1.csv")
# wgs1 = os.path.abspath("ATGC044/Genomes/NC_009881.1.fa")
# wgs2 = os.path.abspath("ATGC044/Genomes/NC_009900.1.fa")
# wgs3 = os.path.abspath("ATGC044/Genomes/NC_017560.1.fa")

#ATGC050 - ready
# csv1 = os.path.abspath("ATGC050/COG_membership/csv/Helicobacter_pylori_83.GCF_000213135.1.csv")
# csv2 = os.path.abspath("ATGC050/COG_membership/csv/Helicobacter_pylori_SouthAfrica7.GCF_000185245.1.csv")
# csv3 = os.path.abspath("ATGC050/COG_membership/csv/Helicobacter_cetorum_MIT_99_5656.GCF_000259275.1.csv")
# wgs1 = os.path.abspath("ATGC050/Genomes/NC_017375.1.fa")
# wgs2 = os.path.abspath("ATGC050/Genomes/NC_017361.1.fa")
# wgs3 = os.path.abspath("ATGC050/Genomes/NC_017735.1.fa")

#ATGC071 - ready
# csv1 = os.path.abspath("ATGC071/COG_membership/csv/Pseudomonas_putida_H8234.GCF_000410575.1.csv")
# csv2 = os.path.abspath("ATGC071/COG_membership/csv/Pseudomonas_putida_NBRC_14164.GCF_000412675.1.csv")
# csv3 = os.path.abspath("ATGC071/COG_membership/csv/Pseudomonas_putida_GB_1.GCF_000019125.1.csv")
# wgs1 = os.path.abspath("ATGC071/Genomes/NC_021491.1.fa")
# wgs2 = os.path.abspath("ATGC071/Genomes/NC_021505.1.fa")
# wgs3 = os.path.abspath("ATGC071/Genomes/NC_010322.1.fa")

#ATGC088 - ready
# csv1 = os.path.abspath("ATGC088/COG_membership/csv/Burkholderia_thailandensis_E264.GCF_000012365.1.csv")
# csv2 = os.path.abspath("ATGC088/COG_membership/csv/Burkholderia_thailandensis_MSMB121.GCF_000385525.1.csv")
# csv3 = os.path.abspath("ATGC088/COG_membership/csv/Burkholderia_pseudomallei_K96243.GCF_000011545.1.csv")
# wgs1 = os.path.abspath("ATGC088/Genomes/NC_007651.1.fa")
# wgs2 = os.path.abspath("ATGC088/Genomes/NC_021173.1.fa")
# wgs3 = os.path.abspath("ATGC088/Genomes/NC_006350.1.fa")

#ATGC089 - ready
# csv1 = os.path.abspath("ATGC089/COG_membership/csv/Burkholderia_cenocepacia_J2315.GCF_000009485.1.csv")
# csv2 = os.path.abspath("ATGC089/COG_membership/csv/Burkholderia_cenocepacia_MC0_3.GCF_000019505.1.csv")
# csv3 = os.path.abspath("ATGC089/COG_membership/csv/Burkholderia_cepacia_GG4.GCF_000292915.1.csv")
# wgs1 = os.path.abspath("ATGC089/Genomes/NC_011000.1.fa")
# wgs2 = os.path.abspath("ATGC089/Genomes/NC_010508.1.fa")
# wgs3 = os.path.abspath("ATGC089/Genomes/NC_018513.1.fa")

#ATGC097 - ready
# csv1 = os.path.abspath("ATGC097/COG_membership/csv/Erwinia_Ejp617.GCF_000165815.1.csv")
# csv2 = os.path.abspath("ATGC097/COG_membership/csv/Erwinia_amylovora_ATCC_49946.GCF_000027205.1.csv")
# csv3 = os.path.abspath("ATGC097/COG_membership/csv/Erwinia_tasmaniensis_Et1_99.GCF_000026185.1.csv")
# wgs1 = os.path.abspath("ATGC097/Genomes/NC_017445.1.fa")
# wgs2 = os.path.abspath("ATGC097/Genomes/NC_013971.1.fa")
# wgs3 = os.path.abspath("ATGC097/Genomes/NC_010694.1.fa")

#ATGC100 - ready
# csv1 = os.path.abspath("ATGC100/COG_membership/csv/Pectobacterium_carotovorum_carotovorum_PC1.GCF_000023605.1.csv")
# csv2 = os.path.abspath("ATGC100/COG_membership/csv/Pectobacterium_carotovorum_carotovorum_PCC21.GCF_000294535.1.csv")
# csv3 = os.path.abspath("ATGC100/COG_membership/csv/Pectobacterium_SCC3193.GCF_000260925.1.csv")
# wgs1 = os.path.abspath("ATGC100/Genomes/NC_012917.1.fa")
# wgs2 = os.path.abspath("ATGC100/Genomes/NC_018525.1.fa")
# wgs3 = os.path.abspath("ATGC100/Genomes/NC_017845.1.fa")

#ATGC104 - ready
# csv1 = os.path.abspath("ATGC104/COG_membership/csv/Bifidobacterium_longum_longum_JCM_1217.GCF_000196555.1.csv")
# csv2 = os.path.abspath("ATGC104/COG_membership/csv/Bifidobacterium_longum_infantis_ATCC_15697_JCM_1222_DSM_20088.GCF_000020425.1.csv")
# csv3 = os.path.abspath("ATGC104/COG_membership/csv/Bifidobacterium_breve_ACS_071_V_Sch8b.GCF_000213865.1.csv")
# wgs1 = os.path.abspath("ATGC104/Genomes/NC_015067.1.fa")
# wgs2 = os.path.abspath("ATGC104/Genomes/NC_011593.1.fa")
# wgs3 = os.path.abspath("ATGC104/Genomes/NC_017218.1.fa")

#ATGC108 - can't find outgroup genome in database ******************

#ATGC111 - ready
# csv1 = os.path.abspath("ATGC111/COG_membership/csv/Aeromonas_hydrophila_hydrophila_ATCC_7966.GCF_000014805.1.csv")
# csv2 = os.path.abspath("ATGC111/COG_membership/csv/Aeromonas_hydrophila_ML09_119.GCF_000401555.1.csv")
# csv3 = os.path.abspath("ATGC111/COG_membership/csv/Aeromonas_salmonicida_salmonicida_A449.GCF_000196395.1.csv")
# wgs1 = os.path.abspath("ATGC111/Genomes/NC_008570.1.fa")
# wgs2 = os.path.abspath("ATGC111/Genomes/NC_021290.1.fa")
# wgs3 = os.path.abspath("ATGC111/Genomes/NC_009348.1.fa")

#ATGC123 - insufficient information to differentiate between species 1 and species 3 ***************

#ATGC125 - missing from database ******************

#ATGC134 - ready
# csv1 = os.path.abspath("ATGC134/COG_membership/csv/Xanthomonas_axonopodis_citri_str_306.GCF_000007165.1.csv")
# csv2 = os.path.abspath("ATGC134/COG_membership/csv/Xanthomonas_axonopodis_citrumelo_F1.GCF_000225915.1.csv")
# csv3 = os.path.abspath("ATGC134/COG_membership/csv/Xanthomonas_oryzae_oryzicola_BLS256.GCF_000168315.3.csv")
# wgs1 = os.path.abspath("ATGC134/Genomes/NC_003919.1.fa")
# wgs2 = os.path.abspath("ATGC134/Genomes/NC_016010.1.fa")
# wgs3 = os.path.abspath("ATGC134/Genomes/NC_017267.2.fa")

#ATGC135 - ready
# csv1 = os.path.abspath("ATGC135/COG_membership/csv/Stenotrophomonas_maltophilia_D457.GCF_000284595.1.csv")
# csv2 = os.path.abspath("ATGC135/COG_membership/csv/Stenotrophomonas_maltophilia_JV3.GCF_000223885.1.csv")
# csv3 = os.path.abspath("ATGC135/COG_membership/csv/Stenotrophomonas_maltophilia_R551_3.GCF_000020665.1.csv")
# wgs1 = os.path.abspath("ATGC135/Genomes/NC_017671.1.fa")
# wgs2 = os.path.abspath("ATGC135/Genomes/NC_015947.1.fa")
# wgs3 = os.path.abspath("ATGC135/Genomes/NC_011071.1.fa")

#ATGC137 - species 1 missing from database ********************

#ATGC138 - ready
# csv1 = os.path.abspath("ATGC138/COG_membership/csv/Francisella_noatunensis_orientalis.GCF_001042565.1.csv")
# csv2 = os.path.abspath("ATGC138/COG_membership/csv/Francisella_philomiragia_philomiragia_ATCC_25017.GCF_000019285.1.csv")
# csv3 = os.path.abspath("ATGC138/COG_membership/csv/Francisella_TX077308.GCF_000219045.1.csv")
# wgs1 = os.path.abspath("ATGC138/Genomes/NZ_CP011923.1.fa")
# wgs2 = os.path.abspath("ATGC138/Genomes/NC_010336.1.fa")
# wgs3 = os.path.abspath("ATGC138/Genomes/NC_015696.1.fa")

#ATGC144 - wrong WGS for top species
# csv1 = os.path.abspath("ATGC144/COG_membership/csv/Borreliella_bissettii_DN127.GCF_000222305.1.csv")
# csv2 = os.path.abspath("ATGC144/COG_membership/csv/Borrelia_burgdorferi_ZS7.GCF_000021405.1.csv")
# csv3 = os.path.abspath("ATGC144/COG_membership/csv/Borrelia_garinii_BgVir.GCF_000239475.1.csv")
# wgs1 = os.path.abspath("ATGC144/Genomes/NC_015907.1.fa")
# wgs2 = os.path.abspath("ATGC144/Genomes/NC_011728.1.fa")
# wgs3 = os.path.abspath("ATGC144/Genomes/NC_017717.1.fa")

#ATGC147 - ready
# csv1 = os.path.abspath("ATGC147/COG_membership/csv/Methanococcus_maripaludis_C5.GCF_000016125.1.csv")
# csv2 = os.path.abspath("ATGC147/COG_membership/csv/Methanococcus_maripaludis_C6.GCF_000018485.1.csv")
# csv3 = os.path.abspath("ATGC147/COG_membership/csv/Methanococcus_maripaludis_C7.GCF_000017225.1.csv")
# wgs1 = os.path.abspath("ATGC147/Genomes/NC_009135.1.fa")
# wgs2 = os.path.abspath("ATGC147/Genomes/NC_009975.1.fa")
# wgs3 = os.path.abspath("ATGC147/Genomes/NC_009637.1.fa")

#ATGC149 - ready
# csv1 = os.path.abspath("ATGC149/COG_membership/csv/Acinetobacter_calcoaceticus_PHEA_2.GCF_000191145.1.csv")
# csv2 = os.path.abspath("ATGC149/COG_membership/csv/Acinetobacter_oleivorans_DR1.GCF_000196795.1.csv")
# csv3 = os.path.abspath("ATGC149/COG_membership/csv/Acinetobacter_baumannii_AB307_0294.GCF_000021145.1.csv")
# wgs1 = os.path.abspath("ATGC149/Genomes/NC_016603.1.fa")
# wgs2 = os.path.abspath("ATGC149/Genomes/NC_014259.1.fa")
# wgs3 = os.path.abspath("ATGC149/Genomes/NC_011595.1.fa")

#ATGC165 - ready
# csv1 = os.path.abspath("ATGC165/COG_membership/csv/Rhodobacter_sphaeroides_ATCC_17029.GCF_000015985.1.csv")
# csv2 = os.path.abspath("ATGC165/COG_membership/csv/Rhodobacter_sphaeroides_KD131.GCF_000021005.1.csv")
# csv3 = os.path.abspath("ATGC165/COG_membership/csv/Rhodobacter_sphaeroides_ATCC_17025.GCF_000016405.1.csv")
# wgs1 = os.path.abspath("ATGC165/Genomes/NC_009049.1.fa")
# wgs2 = os.path.abspath("ATGC165/Genomes/NC_011963.1.fa")
# wgs3 = os.path.abspath("ATGC165/Genomes/NC_009428.1.fa")

#ATGC171 - ready
# csv1 = os.path.abspath("ATGC171/COG_membership/csv/Thermoanaerobacter_X514.GCF_000019065.1.csv")
# csv2 = os.path.abspath("ATGC171/COG_membership/csv/Thermoanaerobacter_wiegelii_Rt8_B1.GCF_000147695.2.csv")
# csv3 = os.path.abspath("ATGC171/COG_membership/csv/Thermoanaerobacter_mathranii_mathranii_str_A3.GCF_000092965.1.csv")
# wgs1 = os.path.abspath("ATGC171/Genomes/NC_010320.1.fa")
# wgs2 = os.path.abspath("ATGC171/Genomes/NC_015958.1.fa")
# wgs3 = os.path.abspath("ATGC171/Genomes/NC_014209.1.fa")

#ATGC177 - ready
# csv1 = os.path.abspath("ATGC177/COG_membership/csv/Prochlorococcus_marinus_AS9601.GCF_000015645.1.csv")
# csv2 = os.path.abspath("ATGC177/COG_membership/csv/Prochlorococcus_marinus_MIT_9301.GCF_000015965.1.csv")
# csv3 = os.path.abspath("ATGC177/COG_membership/csv/Prochlorococcus_marinus_MIT_9312.GCF_000012645.1.csv")
# wgs1 = os.path.abspath("ATGC177/Genomes/NC_008816.1.fa")
# wgs2 = os.path.abspath("ATGC177/Genomes/NC_009091.1.fa")
# wgs3 = os.path.abspath("ATGC177/Genomes/NC_007577.1.fa")

#ATGC181 - ready
# csv1 = os.path.abspath("ATGC181/COG_membership/csv/Caldicellulosiruptor_bescii_DSM_6725.GCF_000022325.1.csv")
# csv2 = os.path.abspath("ATGC181/COG_membership/csv/Caldicellulosiruptor_hydrothermalis_108.GCF_000166355.1.csv")
# csv3 = os.path.abspath("ATGC181/COG_membership/csv/Caldicellulosiruptor_owensensis_OL.GCF_000166335.1.csv")
# wgs1 = os.path.abspath("ATGC181/Genomes/NC_012034.1.fa")
# wgs2 = os.path.abspath("ATGC181/Genomes/NC_014652.1.fa")
# wgs3 = os.path.abspath("ATGC181/Genomes/NC_014657.1.fa")

#ATGC188 - outgroup is missing from database ***********

#ATGC189 - outgroup is missing from database ***********

#ATGC199 - ready
# csv1 = os.path.abspath("ATGC199/COG_membership/csv/Paenibacillus_polymyxa_E681.GCF_000146875.3.csv")
# csv2 = os.path.abspath("ATGC199/COG_membership/csv/Paenibacillus_polymyxa_SC2.GCF_000164985.3.csv")
# csv3 = os.path.abspath("ATGC199/COG_membership/csv/Paenibacillus_terrae_HPL_003.GCF_000235585.1.csv")
# wgs1 = os.path.abspath("ATGC199/Genomes/NC_014483.2.fa")
# wgs2 = os.path.abspath("ATGC199/Genomes/NC_014622.2.fa")
# wgs3 = os.path.abspath("ATGC199/Genomes/NC_016641.1.fa")

#ATGC201 - ready
# csv1 = os.path.abspath("ATGC201/COG_membership/csv/Bartonella_grahamii_as4aup.GCF_000022725.1.csv")
# csv2 = os.path.abspath("ATGC201/COG_membership/csv/Bartonella_tribocorum_CIP_105476.GCF_000196435.1.csv")
# csv3 = os.path.abspath("ATGC201/COG_membership/csv/Bartonella_quintana_Toulouse.GCF_000046685.1.csv")
# wgs1 = os.path.abspath("ATGC201/Genomes/NC_012846.1.fa")
# wgs2 = os.path.abspath("ATGC201/Genomes/NC_010161.1.fa")
# wgs3 = os.path.abspath("ATGC201/Genomes/NC_005955.1.fa")

#ATGC210 - species 2 is missing from database *************

#ATGC213 - ready
# csv1 = os.path.abspath("ATGC213/COG_membership/csv/Methylobacterium_extorquens_AM1.GCF_000022685.1.csv")
# csv2 = os.path.abspath("ATGC213/COG_membership/csv/Methylobacterium_extorquens_PA1.GCF_000018845.1.csv")
# csv3 = os.path.abspath("ATGC213/COG_membership/csv/Methylobacterium_populi_BJ001.GCF_000019945.1.csv")
# wgs1 = os.path.abspath("ATGC213/Genomes/NC_012808.1.fa")
# wgs2 = os.path.abspath("ATGC213/Genomes/NC_010172.1.fa")
# wgs3 = os.path.abspath("ATGC213/Genomes/NC_010725.1.fa")

#ATGC234 - ready
# csv1 = os.path.abspath("ATGC234/COG_membership/csv/Clavibacter_michiganensis_michiganensis_NCPPB_382.GCF_000063485.1.csv")
# csv2 = os.path.abspath("ATGC234/COG_membership/csv/Clavibacter_michiganensis_nebraskensis_NCPPB_2581.GCF_000355695.1.csv")
# csv3 = os.path.abspath("ATGC234/COG_membership/csv/Clavibacter_michiganensis_sepedonicus.GCF_000069225.1.csv")
# wgs1 = os.path.abspath("ATGC234/Genomes/NC_009480.1.fa")
# wgs2 = os.path.abspath("ATGC234/Genomes/NC_020891.1.fa")
# wgs3 = os.path.abspath("ATGC234/Genomes/NC_010407.1.fa")

#ATGC252 - ready
# csv1 = os.path.abspath("ATGC252/COG_membership/csv/Anaeromyxobacter_K.GCF_000020805.1.csv")
# csv2 = os.path.abspath("ATGC252/COG_membership/csv/Anaeromyxobacter_dehalogenans_2CP_1.GCF_000022145.1.csv")
# csv3 = os.path.abspath("ATGC252/COG_membership/csv/Anaeromyxobacter_dehalogenans_2CP_C.GCF_000013385.1.csv")
# wgs1 = os.path.abspath("ATGC252/Genomes/NC_011145.1.fa")
# wgs2 = os.path.abspath("ATGC252/Genomes/NC_011891.1.fa")
# wgs3 = os.path.abspath("ATGC252/Genomes/NC_007760.1.fa")

### FUNCTIONS ###

def main():

    #Get lists of ineligible genes (3' UTR < 30 nts in any of the genomes)
    ineligible1 = check_IGR(csv1)
    ineligible2 = check_IGR(csv2)
    ineligible3 = check_IGR(csv3)

    #Get all of the COG names in all three genomes
    cogs1 = get_cogs(csv1)
    cogs2 = get_cogs(csv2)
    cogs3 = get_cogs(csv3)

    #Count frequency of COGs in each list (only want COGs which appear once per genome)
    fq1 = defaultdict(int)
    for i in cogs1:
        fq1[i] += 1

    fq2 = defaultdict(int)
    for j in cogs2:
        fq2[j] += 1

    fq3 = defaultdict(int)
    for k in cogs3:
        fq3[k] += 1

    #Find which COGs feature once in each dictionary
    useful_cogs = []

    for l in cogs1:
        if l not in ineligible1 and l not in ineligible2 and l not in ineligible3:
            if fq1[l] == 1 and fq2[l] == 1 and fq3[l] == 1:
                useful_cogs.append(l)

    #For each useful COG, extract a fasta sequence for each genome
    for m in useful_cogs:
        section1 = get_fasta(m, csv1, wgs1)
        section2 = get_fasta(m, csv2, wgs2)
        section3 = get_fasta(m, csv3, wgs3)

        section1_seq = section1.split('\n')[1]
        section2_seq = section2.split('\n')[1]
        section3_seq = section3.split('\n')[1]

        if section1_seq != '' and section2_seq != '' and section3_seq != '':
            fasta_total = section1 + section2 + section3

        # Choose a subdirectory, depending on source file used
        subdir = 'Unaligned'
        filename = m + '.txt'
        filepath = os.path.join(subdir, filename)
        f = open(filepath, 'a')
        f.write(fasta_total)
        f.close()


def check_IGR(csv):

    #Set empty list for coords
    info = []
    ineligible = []

    #Open file
    raw_csv = open(csv).read()
    split = raw_csv.split('\n')
    for i in split[1:]:
        if i != '':
            data = i.split(',')
            id = data[1]
            coords = data[5]

            if coords[0:10] == 'complement' and coords[0:15] != 'complement(join':
                c_split = coords.split('..')
                start = int(c_split[0].replace('(', "").replace('complement', '').replace('<', '').replace('>', ''))
                stop = int(c_split[1].replace(')', "").replace('<', '').replace('>', ''))
                complement = '-'
            elif re.findall(r'\d+', coords) == coords.split('..'):
                c_split = coords.split('..')
                start = int(c_split[0].replace('(', "").replace('<', '').replace('>', ''))
                stop = int(c_split[1].replace(')', "").replace('<', '').replace('>', ''))
                complement = '+'
            else:
                continue

            chunk = [id, start, stop, complement]
            info.append(chunk)

    #For each gene, calculate the size of the 3' intergenic space
    for i,n in enumerate(info):
        if info[i][3] == '+':
            if i != len(info)-1:
                IGR = info[i+1][1] - info[i][2]
                if IGR < 30:
                    ineligible.append(info[i][0])

        elif info[i][3] == '-':
            if i != 0:
                IGR = info[i][1] - info[i-1][2]
                if IGR < 30:
                    ineligible.append(info[i][0])

    return ineligible


def get_cogs(csv):

    raw_csv = open(csv).read()
    split = raw_csv.split('\n')
    cog_list = []

    for i in split[1:]:
        if i != '':
            data = i.split(',')
            cog = data[1]
            cog_list.append(cog)

    return cog_list


def get_fasta(cog, csv, wgs):

    #First, get cleaned whole genome seq
    raw_wgs = open(wgs).read()
    split_by_line = raw_wgs.split('\n')
    clean = ''.join(split_by_line[1:])
    fasta_line = split_by_line[0].split(',')[0] + '; ' + cog + ';\n'

    #Next, get coordinates for the appropriate gene
    raw = open(csv).read()
    split = raw.split('\n')

    for i in split:
        if cog in i:
            data = i.split(',')
            coords = data[5].replace('<', '').replace('>', '')

            if coords[0:10] == 'complement' and coords[0:15] != 'complement(join':
                c_split = coords.split('..')
                start = int(c_split[0].replace('(', "").replace('complement', '').replace('<', '').replace('>', ''))
                stop = int(c_split[1].replace(')', "").replace('<', '').replace('>', ''))
                comp = clean[start-19:stop].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()
                reverse_complement = comp[::-1]
                sequence = reverse_complement
            elif re.findall(r'\d+', coords) == coords.split('..'):
                c_split = coords.split('..')
                start = int(c_split[0].replace('(', "").replace('<', '').replace('>', ''))
                stop = int(c_split[1].replace(')', "").replace('<', '').replace('>', ''))
                sequence = clean[start-1:stop+18]
            else:
                sequence = ''

    fasta_section = fasta_line + sequence + '\n'

    return fasta_section

### RUN ###

if __name__ == '__main__':
    main()
