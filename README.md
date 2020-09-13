# CDRF_Metagenomics

[![DOI](https://zenodo.org/badge/126870655.svg)](https://zenodo.org/badge/latestdoi/126870655)

Scripts used for processing of metagenomic data on CDRF project
1. Functional_Annotation --> FOAM_Hmms, ghostKOALA, HMMER and Prokka
2. Identifying_Nitrogen_Fixers --> Custom python scripts to extract information out of tables from database.
3. Sequence_QC_KneadData --> Shell scripts to run kneaddata for QC of metagenomic reads.
4. Taxonomy_Kaiju --> Shell scripts to run Kaiju for taxanomic assignments.
5. Anvio_Workflow  --> Scripts for constructing anvio database of contigs and their taxonomy and functional annotation. 
6. Assembly --> Shell scripts for digital normalization and assembly of reads.
7. Assembly by Farm  --> Assembly by groups.
8. Nitro_Hmmscan_to_GeneCount.py is a python script that loops through tblout output files from hmmscan (of nitrogen genes) in a directory and converts an output file of counts of all genes with a high at a certain E-value. 
9. Resfam_Hmmscan_to_GeneCount.py is a python script that loops through tblout output files from hmmscan (Resfam databasae) in a directory and converts an output file of counts of all genes with a high at a certain E-value. 
