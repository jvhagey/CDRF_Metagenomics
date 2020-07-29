#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019

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
    parser = ArgumentParser(prog="anvi-get-taxonomy-for-genecall.py", description="Parses anvi-export-table --table files")
    parser.add_argument("-t", "--taxa-level",type=str, dest="taxalevel", action="store", default=None, help='Pick taxa level you want reported')
    parser.add_argument("--print-gene", dest='printgene',action="store_true", default=False, help="Print contigs with your choose combination of taxa and gene call(s).")
    parser.add_argument("-i", "--indir", dest="indirname",action="store", default=None, required=True, help="Input directory name were taxon_names, genes_taxonomy and hmm_hits .txt files are found (required)")
    parser.add_argument("-o", "-outfile", dest="outfile", action="store", default=None, help="Name of file to write output too.")
    parser.add_argument("-g", "--gene-call", type=str, dest="gene", action="store", default=None, required=False, help="Name of gene to search in dataframe for. You can pass multiple genes by separating them by a pipe. Ex: 'geneA|geneB'")
    args = parser.parse_args()
    return args

#set colors for warnings so they are noticable
CRED = '\033[91m'+'\nWarning:'
CYEL = '\033[93m'+'\nJust FYI:'
CPUR = '\033[95m'+'\nQuick Question:'
CEND = '\033[0m'

def check_gene(gene,dataframe):
    genes = dataframe['gene_callers_id'].astype(str).unique().tolist()
    if gene == "list" or None:
        print(CRED +"\nHere is list to pick from:")
    elif gene is not None and not any(re.findall(r'|'.join(dataframe['gene_callers_id'].astype(str)),gene, re.IGNORECASE)):
        print(CRED + "\nIt looks like this gene isn't in the dataframe. :(" )
        Q = input(CPUR +'\nYikes there are ' + genes.value_counts() + ' unique genes here do you want them all printed out (yes/no)?' + CEND).lower()
        if Q.startswith('y'):#check to see if this is the dataframe you want printed
            print(CRED +"\nHere is list to pick from:")
            print(genes)
        if Q.startswith('n'):
            exit()

def main ():
    args = parse_cmdline() #Parse command-line to get args
    #Calling in files
    os.chdir(args.indirname) #change directory
    Anvi_table_parser_new3.check_taxa_input(args.taxalevel) #check to make sure taxonmy argument is valid
    #check that files are present
    Anvi_table_parser_new3.check_files_existance(args.indirname)
    #run_all_works2.get_datatables(args.indirname)
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
    check_gene(args.gene, new3)
    df = new3[new3['gene_callers_id'].isin([args.gene])]
    print(CYEL + "\n")
    print(df['t_'+args.taxalevel].unique().tolist()) #get list of unique taxa that contain nifH.
    Q = input(CPUR +'\nDo you want to print this data frame printed to a csv file (yes/no)?' + CEND).lower()
    if Q.startswith('y'):#check to see if this is the dataframe you want printed
        if outfile is None: #check that outfile name was given
            print(CRED +"\n You need to select an out file name to write dataframe to!\n"+ CEND)
            exit()
        if outfile is not None:
            with open(outfile, "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(data_frame.to_csv(outfile, header=0, sep='\t', mode='w', index=None))
                print(CYEL + "\nCool we just printed your data frame to the file "+outfile +CEND)
                exit()

#['gene_callers_id', 'taxon_id', 't_phylum', 't_class', 't_order', 't_family', 't_genus', 't_species', 'entry_id', 'source_x', 'gene_name', 'gene_hmm_id', 'e_value', 'contig', 'partial']

if __name__ == '__main__':
    main()
