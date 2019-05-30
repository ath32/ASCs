# ASCs

================

This repository contains the source code for running the analysis described in the paper. Not all of the scripts in this repository will be relevant to the complete replication of the study - if you wish to do this, please send me an e-mail and we can sort something out (a.t.ho@bath.ac.uk).

###Instructions:
- All files run from the root directory.
- Python scripts are found in '0_Python_scripts' (all scripts are run using Python3)
- R scripts are found in '0_R_scripts'
- CSVs produced by Python scripts can be found in 'Author CSVs'

###The folders found in my working directory are labelled 0-4, all available on request:
- 0 - Scripts
- 1 - EMBL files 
- 2 - FASTA files 
- 3 - PaxDB raw data and appropriate EMBL genomes
- 4 - Outputs (includes CSVs)

###Python scripts:

Python scripts are categorised 1-11 as follows. A spreadsheet with more detailed script information (Python and R) can be found in the document 'Scripts_Explained'.

- 1 - Filtering and obtaining FASTA files
- 2 - Dinucleotide-controlled simulation experiment
- 3 - Obtain ASC frequency information - includes relative usage
- 4 - Mollicutes analysis
- 5 - Subset analysis e.g. TGA HEGs vs TAA LEGs
- 6 - Third stops analysis
- 7 - +4T identification and tests
- 8 - Primary stop usage
- 9 - +4T-controlled simulation experiment
- 10 - T-codon enrichment scores
- 11 - Lengths and expression association

*** Codon switch analysis files (including scripts) are found in the 'COGs' folder ***

###FASTA folders/files explained:

My working repository contains many different folders of FASTA files (I am happy to send these to you upon request). These are each explained below:

2_FASTA_CDS_all: This folder contains FASTA files for our filtered bacteria dataset. Sequences contained in the FASTA files are coding sequence information for all qualifying genes that satisfy our filters.

2_FASTA_CDS_expression: This folder contains FASTA files for all genomes for which protein expression data was available. Sequences contained in the FASTA files are coding sequence information for all qualifying genes that satisfy our filters. Files in this folder contain protein expression data for each gene.

2_FASTA_CDS_hegs: This folder contains FASTA files for all genomes for which protein expression data was available. Sequences contained in the FASTA files are coding sequence information for all highly expressed genes (see methods for how these are defined).

2_FASTA_CDS_legs: This folder contains FASTA files for all genomes for which protein expression data was available. Sequences contained in the FASTA files are coding sequence information for all lowly expressed genes (see methods for how these are defined).

2_FASTA_Eubacteria_cds_mollicutes: This folder contains FASTA files for all mollicute genomes found in the unfiltered bacterial genomes dataset. Sequences contained in the FASTA files are the UTRs for all qualifying genes.

2_FASTA_Eubacteria_cds_TT4: This folder contains FASTA files for translation table 4 genomes from our filtered bacteria dataset. Sequences contained in the FASTA files are the UTRs for all qualifying genes.

2_FASTA_Eubacteria_cds_TT11: This folder contains FASTA files for translation table 11 genomes from our filtered bacteria dataset. Sequences contained in the FASTA files are the UTRs for all qualifying genes.

2_FASTA_Eubacteria_cds_TT11_TAA: This folder contains FASTA files for translation table 11 genomes from our filtered bacteria dataset. Sequences contained in the FASTA files are the UTRs for all qualifying genes that terminate TAA.

2_FASTA_Eubacteria_cds_TT11_TAG: This folder contains FASTA files for translation table 11 genomes from our filtered bacteria dataset. Sequences contained in the FASTA files are the UTRs for all qualifying genes that terminate TAG.

2_FASTA_Eubacteria_cds_TT11_TGA: This folder contains FASTA files for translation table 11 genomes from our filtered bacteria dataset. Sequences contained in the FASTA files are the UTRs for all qualifying genes that terminate TGA.

2_FASTA_PaxDb: This folder contains FASTA files for all genomes for which protein expression data was available. Sequences contained in the FASTA files are the UTRs for all qualifying genes that satisfy our filters.

2_HEGs_fasta: This folder contains FASTA files for all genomes for which protein expression data was available. Sequences contained in the FASTA files are the UTRs for all qualifying highly expressed genes (see methods for how these are defined).

2_LEGs_fasta: This folder contains FASTA files for all genomes for which protein expression data was available. Sequences contained in the FASTA files are the UTRs for all qualifying lowly expressed genes (see methods for how these are defined).


### Alexander Ho
### Milner Centre for Evolution
### University of Bath
### a.t.ho@bath.ac.uk
