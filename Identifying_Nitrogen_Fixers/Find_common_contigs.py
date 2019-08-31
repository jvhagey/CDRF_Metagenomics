#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019
## Given files extracted from an anvio.db and Anvi_table_parser.py script will output common taxa or contigs found in common between multiple methods

import Anvi_table_parser_new3
from subprocess import Popen,PIPE
import pandas as pd
import os, sys, re
from argparse import ArgumentParser

def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="Identify_Nif_Taxonomy.py", description="Parses anvi-export-table --table files")
    parser.add_argument("-i", "--indir", dest="indirname",action="store", default=None, required=True, help="Input directory name were taxon_names, genes_taxonomy and hmm_hits .txt files are found (required)")
    parser.add_argument("-t", "--taxa-level",type=str, dest="taxalevel", action="store", default=None, help='Pick taxa level you want reported')
    parser.add_argument("--print-taxa", dest='printtaxa',action="store_true", default=False, help="Print taxa with nifHDK genes at level chosen to text file.")
    parser.add_argument("--get-taxa-contigs", dest='gettaxacontigs',action="store_true", default=False, help="Get contigs with nifHDK genes to text file.")
    parser.add_argument("-g", "--gene", type=str, dest="gene", action="store", default=None, required=False, help="Name of gene to search in dataframe for. You can pass multiple genes by separating them by a pipe. Ex: 'geneA|geneB'")
    parser.add_argument("--get-gene-contigs", dest='getgenecontigs',action="store_true", default=False, help="Get contigs containing gene(s) of choice to text file.")
    args = parser.parse_args()
    return args

#set colors for warnings so they are noticable
CRED = '\033[91m'+'\nWarning:'
CYEL = '\033[93m'
CPUR = '\033[95m'+'\nQuick Question:'
CEND = '\033[0m'

def arg_check(printtaxa, gettaxacontigs, gene):
    if printtaxa == True and gettaxacontigs == True:
        print(CRED + "\n You can't pass both --get-contigs and --print-taxa. Choose just one at a time :(")
        exit()
    if printtaxa == False and gettaxacontigs == False and gene is None:
        print(CRED + "\n You need to pick a mode, either pass gene name with '-g gene_x' or --get-contigs or --print-taxa. Choose just one at a time :(")
        exit()

def get_datatables(indirname):
        Anvi_table_parser_new3.check_files_existance(indirname)
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
        #"Nif_Hmms"
        #calling in dataframe that has gene names for Srijak's custom nif gene hmms.
        df_srijak_hmm_names=pd.read_csv("Srijak_hmm_names.txt", sep='\t', header=0)
        #searching for particular rows with matching source_name
        global df_nif_hmms2, df_foam_hmms4, GK, Cog
        df_nif_hmms = new3[new3['source_x'].str.contains("Nif_Hmms")]
        df_nif_hmms2 = pd.merge(df_nif_hmms, df_srijak_hmm_names, on='gene_name')
        #"Foam_Nitro":
        #loading dataframes
        KO_modules = pd.read_csv("KO_modules.txt", sep='\t', header=0) #change column name to allow merging later
        df_Foam_names = pd.read_csv("FOAM-onto_rel1.txt", sep='\t', header=0) #table with KO number and FOAM level assignment
        #searching for particular rows with matching source_name
        df_foam_hmms = new3[new3['source_x'].str.contains("Foam")]
        df_foam_hmms['KO'] = df_foam_hmms.gene_name.str.extract('KO:([^_]*)', expand=False) #capture everything between KO: and _
        df_foam_hmms2 = pd.merge(df_foam_hmms, df_Foam_names, on='KO') #add foam level assignment info to contigs with hits from foam hmms
        df_foam_hmms3 = df_foam_hmms2[df_foam_hmms2['L1'].str.contains("11_Nitrogen cycle")]
        df_foam_hmms4 = pd.merge(df_foam_hmms3, KO_modules, on='KO')
        # "KeggGhostKoala"
        df_kegg_hmm_names = pd.read_csv("KeggOrthology_Nitrogen_Table1.txt", sep=',', header=0) #table with KeggOrthology info
        df_kegg_hmm_names = df_kegg_hmm_names.rename(columns={'accession':'KO'})
        df_genfun = pd.read_csv("gene_functions.txt", sep='\t', header=0) #has gene_callers_id and the hmm hit info
        df_GK = df_genfun[df_genfun['source'].str.contains("KeggGhostKoala")]
        df_GK_gen = pd.merge(df_GK, df_genes_contigs, on='gene_callers_id') #combine to get contig information with hmm hit info
        GK = pd.merge(new, df_GK_gen, on='gene_callers_id')
        GK = GK.drop(columns=['start','stop','direction','taxon_id','source_y','version', 'entry_id', 'partial'])
        #"COG"
        df_Cog = df_genfun[df_genfun['source'].str.contains("COG")]
        df_Cog_gen = pd.merge(df_Cog, df_genes_contigs, on='gene_callers_id') #combine to get contig information with hmm hit info
        Cog = pd.merge(new, df_Cog_gen, on='gene_callers_id')
        Cog = Cog.drop(columns=['start','stop','direction','taxon_id','source_y','version', 'entry_id', 'partial'])
def get_gene_df(gene, getgenecontigs):
    global Foam_gene_list, GK_gene_list, Cog_gene_list, Nif_Hmms_gene_list
    if gene is not None:
        Foam_gene = df_foam_hmms4[df_foam_hmms4['Gene'].str.contains(gene, flags=re.IGNORECASE)]
        Foam_gene_list = Foam_gene.contig.unique().tolist()
        #print(CYEL + "\n This is a list of contigs found with Foam Hmms for "+gene+":\n")
        #print(Foam_gene.contig.unique().tolist())
        #global GK_gene_list
        GK_gene = GK[GK['function'].str.contains(gene, flags=re.IGNORECASE)]
        GK_gene_list = GK_gene.contig.unique().tolist()
        #print(CYEL + "\n This is a list of contigs found with GhostKOALA for "+gene+":\n")
        #print(GK_gene.contig.unique().tolist())
        #global Cog_gene_list
        Cog_gene = Cog[Cog['function'].str.contains(gene, flags=re.IGNORECASE)]
        Cog_gene_list = Cog_gene.contig.unique().tolist()
        #print(CYEL + "\n This is a list of contigs found with Cogs for "+gene+":\n")
        #print(Cog_gene.contig.unique().tolist())
        #global Nif_Hmms_gene_list
        Nif_Hmms_gene = df_nif_hmms2[df_nif_hmms2['Gene_description'].str.contains(gene, flags=re.IGNORECASE)]
        Nif_Hmms_gene_list = Nif_Hmms_gene.contig.unique().tolist()
        #print(CYEL + "\n This is a list of contigs found with nif TigrFams for "+gene+":\n")
        #print(Nif_Hmms_gene.contig.unique().tolist())
        if getgenecontigs == True:
            with open('contigs_Nif_Hmms_'+gene+'.txt', "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(Nif_Hmms_gene_list.to_csv('contigs_Nif_Hmms_'+gene+'.txt', header=True, sep='\t', mode='w', index=None))
            with open('contigs_Cog_'+gene+'.txt', "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(Cog_gene_list.to_csv('contigs_Cog_'+gene+'.txt', header=True, sep='\t', mode='w', index=None))
            with open('contigs_GK_'+gene+'.txt', "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(GK_gene_list.to_csv('contigs_GK_'+gene+'.txt', header=True, sep='\t', mode='w', index=None))
            with open('contigs_Foam_Hmms_'+gene+'.txt', "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(Foam_gene_list.to_csv('contigs_Foam_Hmms_'+gene+'.txt', header=True, sep='\t', mode='w', index=None))
            print(CYEL + "\nCool we just printed some contigs into individual files. Now you can search for these contigs in MAG bins!")
            exit()
        else:
            pass
def get_taxa(taxalevel, printtaxa, gettaxacontigs):
    global Foam_taxa_list, GK_taxa_list, Cog_taxa_list, Nif_Hmms_taxa_list
    if printtaxa == True:
        Foam_taxa_list = Anvi_table_parser_new3.find_HDK(df_foam_hmms4, "Foam_Nitro", True, None, taxalevel, 'nifHDK_Foam_Nitro_'+taxalevel+".txt", False, printtaxa, gettaxacontigs)
        #global GK_taxa_list
        GK_taxa_list = Anvi_table_parser_new3.find_HDK(GK, "KeggGhostKoala", True, None, taxalevel, 'nifHDK_GK_'+taxalevel+".txt", False, printtaxa, gettaxacontigs)
        #global Cog_taxa_list
        Cog_taxa_list = Anvi_table_parser_new3.find_HDK(Cog, "COG", True, None, taxalevel, 'nifHDK_Cog_'+taxalevel+".txt", False, printtaxa, gettaxacontigs)
        #global Nif_Hmms_taxa_list
        Nif_Hmms_taxa_list = Anvi_table_parser_new3.find_HDK(df_nif_hmms2, "Nif_Hmms", True, None, taxalevel, 'nifHDK_nif_hmms_'+taxalevel+".txt", False, printtaxa, gettaxacontigs)
    else:
        Foam_taxa_list = Anvi_table_parser_new3.find_HDK(df_foam_hmms4, "Foam_Nitro", True, None, taxalevel, None, False, printtaxa, gettaxacontigs)
        #global GK_taxa_list
        GK_taxa_list = Anvi_table_parser_new3.find_HDK(GK, "KeggGhostKoala", True, None, taxalevel, None, False, printtaxa, gettaxacontigs)
        #global Cog_taxa_list
        Cog_taxa_list = Anvi_table_parser_new3.find_HDK(Cog, "COG", True, None, taxalevel, None, False, printtaxa, gettaxacontigs)
        #global Nif_Hmms_taxa_list
        Nif_Hmms_taxa_list = Anvi_table_parser_new3.find_HDK(df_nif_hmms2, "Nif_Hmms", True, taxalevel, None, False, printtaxa, gettaxacontigs)
def find_common(nif_hmms, foam_hmms, GK, cog, printtaxa, gettaxacontigs, gene):
    if printtaxa == True:
        word = "taxa"
    if gettaxacontigs == True:
        word = "contigs"
    if gene is not None:
        word = "contigs for " + gene
    # Converting the arrays into sets
    s1 = set(nif_hmms)
    s2 = set(foam_hmms)
    s3 = set(GK)
    s4 = set(cog)
    # Calculates intersection of sets
    set1 = s1.intersection(s2)
    print(CYEL + "\nThe " +word+ " that were found in common between custom TigrFams and foam_hmms:")
    print(list(set1))
    set2 = s3.intersection(s4)
    print("\nThe " +word+ " that were found in common between GK and Cogs:")
    print(list(set2))
    final_list = set1.intersection(set2)
    print("\nThe " +word+ " that were found in common will all three methods are:")
    print(list(final_list))

def main():
    args = parse_cmdline()
    os.chdir(args.indirname) #change directory
    arg_check(args.printtaxa, args.gettaxacontigs, args.gene)
    get_datatables(args.indirname)
    Anvi_table_parser_new3.check_gene_input(df_foam_hmms4, "Foam_Nitro", args.gene)
    Anvi_table_parser_new3.check_gene_input(df_nif_hmms2, "Nif_Hmms", args.gene)
    get_gene_df(args.gene, args.getgenecontigs)
    if args.gene is not None:
        find_common(Nif_Hmms_gene_list, Foam_gene_list, GK_gene_list, Cog_gene_list, args.printtaxa, args.gettaxacontigs, args.gene)
    if args.printtaxa == True or args.gettaxacontigs == True:
        get_taxa(args.taxalevel, args.printtaxa, args.gettaxacontigs)
        find_common(Nif_Hmms_taxa_list,Foam_taxa_list,GK_taxa_list,Cog_taxa_list,args.printtaxa, args.gettaxacontigs, None)

    #print(list(set(Foam_taxa_list).intersection(GK_taxa_list)))

if __name__ == '__main__':
    main()
