#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2020
## Given files extracted from an anvio.db script will return the number of gene calls in gene

#importing packages
import glob, os, re, csv
import pandas as pd
import re
import csv
import Anvi_table_parser_new3
from functools import reduce
from argparse import ArgumentParser
#import feather

#Defining function to process command-line arguments
def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="get_gene_calls_from_contig.py", description="Parses anvi-export-table --table files")
	parser.add_argument("-i", "--indir", dest="indirname",action="store", default=None, required=True, help="Input directory name were taxon_names, genes_taxonomy and hmm_hits .txt files are found (required)")
	parser.add_argument("-o", "-outfile", dest="outfile", action="store", default=None, help="Name of file to write output too.")
	parser.add_argument("-c", "--contig-ID", type=str, dest="contig", action="store", default=None, required=False, help="Contig ID you want to search for in dataframe. You can pass multiple contigs by separating them by a pipe. Ex: 'contig#|contig#'")
	args = parser.parse_args()
	return args

#set colors for warnings so they are noticable
CRED = '\033[91m'+'\nWarning:'
CYEL = '\033[93m'+'\nJust FYI:'
CPUR = '\033[95m'+'\nQuick Question:'
CEND = '\033[0m'

def main ():
	args = parse_cmdline() #Parse command-line to get args
	#Calling in files
	os.chdir(args.indirname) #change directory
	Anvi_table_parser_new3.check_taxa_input(args.taxalevel) #check to make sure taxonmy argument is valid
	#check that files are present
	Anvi_table_parser_new3.check_files_existance(args.indirname)
	#making dataframes
	df_genes_contigs = pd.read_csv("genes_in_contigs.txt", sep='\t', header=0) #has gene_callers_id and contig information
	df_hmm_hits = pd.read_csv("hmm_hits.txt", sep='\t', header=0) #has gene_callers_id and the hmm hit info
	df_gene_taxa = pd.read_csv("genes_taxonomy.txt", sep='\t', header=0) #has gene_callers_id and taxon_id
	df_taxa_names = pd.read_csv("taxon_names.txt", sep='\t', header=0)  #has taxon_id and their assigned taxa
	#Merging dataframes
	new = pd.merge(df_gene_taxa, df_taxa_names, on='taxon_id') #combine to get gene_callers_id to assigned taxa
	new2 = pd.merge(df_hmm_hits, df_genes_contigs, on='gene_callers_id') #combine to get contig information with hmm hit info
	new3 = pd.merge(new, new2, on='gene_callers_id') #combine to get taxanomic assignment of contigs and hmm hits in those contigs
	new3 = new3.drop(columns=['start','stop','direction','gene_unique_identifier','source_y','version']) #droping columns
	df = new3[new3['contig'].isin([args.contig])]
	print(CYEL +"\n Here are the gene calls in this contig!\n"+ CEND)
	print(df['gene_callers_id'].unique().tolist()) #get list of unique gene calls in a contig.
	print(CYEL +"\n This is the number of gene calls in this contig!\n"+ CEND)
	print(len(df['gene_callers_id'].unique().tolist()))
	Q = input(CPUR +'\nDo you want to print this data frame printed to a csv file (yes/no)?' + CEND).lower()
	if Q.startswith('y'):#check to see if this is the dataframe you want printed
		if outfile is None: #check that outfile name was given
			print(CRED +"\n You need to select an out file name to write dataframe to!\n"+ CEND)
			exit()
		if outfile is not None:
			with open(outfile, "w") as text_file: #writing out list of gene calls in contig
				pd.DataFrame(data_frame.to_csv(outfile, header=0, sep='\t', mode='w', index=None))
				print(CYEL + "\nCool we just printed your data frame to the file "+outfile +CEND)
				exit()

#['gene_callers_id', 'taxon_id', 't_phylum', 't_class', 't_order', 't_family', 't_genus', 't_species', 'entry_id', 'source_x', 'gene_name', 'gene_hmm_id', 'e_value', 'contig', 'partial']

if __name__ == '__main__':
	main()
