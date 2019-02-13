# ASCs

================

[Pre-amble] This repository contains the source code for running the analysis described in the paper.

Instructions:
- All files run from the root directory.
- Scripts 1,2 are used to download the EMBL files. [Genome download information]
- A list of accessions for the genomes used can be found in accession_lists.
- Eukaryote genomes were downloaded directly from Ensembl and need to be placed in a root folder '1_Eubacteria_filtered_genomes'
- The following empty folders should be created in the root directory: [Enter in here]
- All python scripts are run using Python 3.

Scripts:


1_bacterial_genome_download
get_genomes.tcl: download genomes (command requires path to accession list)
group_genomes.py: group genomes into one folder
2_sort_genomes
_sortFiles.py: sort the genomes
3_parse_files
_parseFiles.py: parse the embl files and output to fasta like formats
_filterGenes.py: filter genes
_parseFiles_t4.py: parse the files for all table 4 genomes
4_ratio_testing
_calculate_ratios.py: calculate enrichment ratios for given site
_chitest.py: chi square test for given site
_site_4_no_overlap.py: calculate enrichment ratios when cds with an overlap to the next cds are discounted
_chitest_site_4_no_overlap: chitest for site 4 with no overlaps
_a_usage_start_cds.py: get the A usage at sites in the 5' domain
_all_t4_genomes.py: calculate fourth site ratios for all table 4 genomes
5_gc_conservation
_codon_nucleotide_proportions.py: calculate the nucleotide proportions at each position
_codon_gc_proportions.py: calculate the gc content at each position
_codon_gc_varaiance.py: calcaulte the gc variance at each position
6_aod_calculation
_aod_second_amino.py: calculate aod scores for amino acids
7_stop_codons
_stopDistanceGenome.py: distances to first +1 stop
_stopDistanceGenomeSecondStop.py: distances to second +1 stop
_stopDistanceGenomeThirdStop.py: distances to third +1 stop
8_expression
_getHighExpGenes.py: get genomes with annotations of the highly expressed genes
_genomeGenes.py: create files with all genes, and the highly expressed genes
_runCodonW.py: calculate CAI scores
_CAI_analysis.py: run analysis on CAI values
_CAI_stats.py: get CAI for the different start codons
9_related_species
_compare_base_changes_related_species.py: compare orthologs between e. coli and shigella
10_leader_genes
_leaderGenes.py: search for upstream leader genes
_filteredLeader.py: filter the genes that may have a leader
_leaderAnalysis.py: analyse the leaders
_leader_leaderless_a_content: get the a content dependent on leader status
11_upstream_cds
_upstream_cds.py: see whether an upstream cds affects fourth site content
12_met_pairs
_met_pairs.py: compare use of amino acids after methionine
13_sd_sequences:
_anti_sd_sequences.py: get the 16s anti sd sequence
_sdSequences.py: calculate the binding strength of cds with the antisd sequence
_analyse_sd_sequences.py: analysis of sd sequences
14_archaea
14.1_get_genomes.tcl: download archaea genomes (command requires path to accession list)
14.2_group_genomes.py: group archaea genomes
14.3_sortFiles.py: sort files
14.4_parseFiles.py: extract to fasta like format
14.5_filterGenes.py: filter the cdss
14.6_calculate_ratios: calculate enrichment ratios for given site
14.7_chitest.py: chitest for site
15_eukaryotes
15.1_filterGenes.py: filter the cdss
15.2_calculate_ratios.py: calculate ratios for given site
16_second_amino_usage
_second_amino_relative_usage.py: amino acid use in the second position
_second_amino.py: get enrichment ratios of amino acids
17_multivar
17_multivar.py: sort outputs for multivar analysis
18_protists
get_genomes.py: download protist genomes
protists.py: calculate enrichment ratios for protists
table_4_tests
Two independent scripts verifying fourth site calculations

