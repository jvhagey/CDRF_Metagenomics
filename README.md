# CDRF_Metagenomics
Scripts used for processing of metagenomic data on CDRF project
1. Sequence_QC_KneadData --> Shell scripts to run kneaddata for QC of metagenomic reads 
2. Taxonomy_Kaiju --> Shell scripts to run Kaiju for taxanomic assignments
3. Assembly --> Shell scripts for digital normalization and assembly of reads
4. Nitro_Hmmscan_to_GeneCount.py is a python script that loops through tblout output files from hmmscan (of nitrogen genes) in a directory and converts an output file of counts of all genes with a high at a certain E-value. 
5. Resfam_Hmmscan_to_GeneCount.py is a python script that loops through tblout output files from hmmscan (Resfam databasae) in a directory and converts an output file of counts of all genes with a high at a certain E-value. 
