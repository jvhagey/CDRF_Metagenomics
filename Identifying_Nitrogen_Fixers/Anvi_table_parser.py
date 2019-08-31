#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019
## Given files extracted from an anvio.db, script will give a .feather file of KO counts to read into R

#importing packages
import glob, os, re, csv
import pandas as pd
import re
import csv
#import feather
from functools import reduce
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="Identify_Nif_Taxonomy.py", description="Parses anvi-export-table --table files")
    parser.add_argument("-t", "--taxa-level",type=str, dest="taxalevel", action="store", default=None, help='Pick taxa level you want reported')
    parser.add_argument("-l", "--foam-level",dest="Foamlevel", action="store", default=None, help="For list of foam level hits use '-l list' as arguement default none.")
    parser.add_argument("--print-contigs", dest='printcontigs',action="store_true", default=False, help="Print contigs with your choose combination of source, taxa and gene.")
    #parser.add_argument("--make-datatables", dest='datatables',action="store_true", default=False, help="Supress questions and don't exit. Only use when importing module.")
    parser.add_argument("--find-nifHDK-in-contigs", dest='nifgenes',action="store_true", default=False, help="Print contigs with nifHDKENB genes.")
    parser.add_argument("-st", "--specific-taxa", dest="specifictaxa", action="store", default=None, help="Give name of specific taxa you are interested in. You need to pass -t for the appropriate level as well.")
    parser.add_argument("-i", "--indir", dest="indirname",action="store", default=None, required=True, help="Input directory name were taxon_names, genes_taxonomy and hmm_hits .txt files are found (required)")
    parser.add_argument("-o", "-outfile", dest="outfile", action="store", default=None, help="Name of file to write output too.")
    parser.add_argument("-g", "--gene", type=str, dest="gene", action="store", default=None, required=False, help="Name of gene to search in dataframe for. You can pass multiple genes by separating them by a pipe. Ex: 'geneA|geneB'")
    parser.add_argument("-s", "--source", type=str, dest="source", action="store", default=None, required=True, help="Indicate source, for list of sources use '-s list' as arguement default none.")
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
    #check that files are present
    check_files_existance(args.indirname)
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
    check_source_input(new3, args.source) #Checking to make sure arguments are one of the options
    check_taxa_input(args.taxalevel) #check to makre sure taxonmy argument is valid
    if args.source == "Nif_Hmms":
        #calling in dataframe that has gene names for Srijak's custom nif gene hmms.
        df_srijak_hmm_names=pd.read_csv("Srijak_hmm_names.txt", sep='\t', header=0)
        #searching for particular rows with matching source_name
        df_nif_hmms = new3[new3['source_x'].str.contains("Nif_Hmms")]
        df_nif_hmms2 = pd.merge(df_nif_hmms, df_srijak_hmm_names, on='gene_name')
        check_gene_input(df_nif_hmms2, args.source, args.gene)
        find_HDK(df_nif_hmms2, args.source, args.nifgenes, args.gene, args.taxalevel, args.outfile, args.printcontigs, False, False)
        if args.gene is None:
            specific_taxa(df_nif_hmms2, args.specifictaxa, args.taxalevel)
            print(CYEL + "\nThis is the full dataframe, of contigs that had hits with from hmms from nifHDKENB etc. \nTo reduce the size of the dataframe pick a gene with or without a taxalevel to look at.\n" + CEND)
            print_exit(df_nif_hmms2, args.printcontigs, args.outfile, args.taxalevel)
            print_taxa(df_nif_hmms2, args.taxalevel, args.outfile)
        if args.gene is not None:
            df_gene_srijak = df_nif_hmms2[df_nif_hmms2['Gene_description'].str.contains(args.gene)] #searching for particular rows with matching gene_name
            specific_taxa(df_gene_srijak, args.specifictaxa, args.taxalevel)
            print(CYEL + "\nThis is a dataframe, of contigs that had hits with from the hmm(s) for your choosen gene(s). \nTo find out what taxa the contigs belong to pass -t with a taxalevel. If you already passed -t just hang on and we will get there\n" + CEND)
            print_exit(df_gene_srijak, args.printcontigs, args.outfile, args.taxalevel)
            print_taxa(df_gene_srijak, args.taxalevel, args.outfile)
    if args.source == "Foam_Nitro":
        #loading dataframes
        KO_modules = pd.read_csv("KO_modules.txt", sep='\t', header=0) #change column name to allow merging later
        df_Foam_names = pd.read_csv("FOAM-onto_rel1.txt", sep='\t', header=0) #table with KO number and FOAM level assignment
        #searching for particular rows with matching source_name
        df_foam_hmms = new3[new3['source_x'].str.contains("Foam")]
        df_foam_hmms['KO'] = df_foam_hmms.gene_name.str.extract('KO:([^_]*)', expand=False) #capture everything between KO: and _
        df_foam_hmms2 = pd.merge(df_foam_hmms, df_Foam_names, on='KO') #add foam level assignment info to contigs with hits from foam hmms
        df_foam_hmms3 = df_foam_hmms2[df_foam_hmms2['L1'].str.contains("11_Nitrogen cycle")]
        df_foam_hmms4 = pd.merge(df_foam_hmms3, KO_modules, on='KO')
        #df_foam_hmms4 = pd.merge(df_foam_hmms3, df_kegg_hmm_names, on='KO') #add in KeggOrthology info, when you do this the dataframe is expanded due to KOs being in more than one category
        check_foamlevel_input(args.Foamlevel)
        check_gene_input(df_foam_hmms4,args.source, args.gene)
        find_HDK(df_foam_hmms4, args.source, args.nifgenes, args.gene, args.taxalevel, args.outfile, args.printcontigs, False, False)
        if args.Foamlevel is None:
            run_gene(df_foam_hmms4, args.gene, args.source, args.specifictaxa, args.taxalevel, args.printcontigs, args.outfile)
        if args.Foamlevel is not None:
            df_foam_level = df_foam_hmms4[df_foam_hmms4['L2'].str.contains(args.Foamlevel)]
            run_gene(df_foam_level, args.gene, args.source, args.specifictaxa, args.taxalevel, args.printcontigs, args.outfile)
        clean_up("df_foam_hmms4.txt")
    if args.source == "KeggGhostKoala" or args.source == "COG":
        #loading dataframes
        df_kegg_hmm_names = pd.read_csv("KeggOrthology_Nitrogen_Table1.txt", sep=',', header=0) #table with KeggOrthology info
        df_kegg_hmm_names = df_kegg_hmm_names.rename(columns={'accession':'KO'})
        df_genfun = pd.read_csv("gene_functions.txt", sep='\t', header=0) #has gene_callers_id and the hmm hit info
        df_Cog_GK = df_genfun[df_genfun['source'].str.contains(args.source)]
        df_Cog_GK_gen = pd.merge(df_Cog_GK, df_genes_contigs, on='gene_callers_id') #combine to get contig information with hmm hit info
        Cog_GK = pd.merge(new, df_Cog_GK_gen, on='gene_callers_id')
        Cog_GK = Cog_GK.drop(columns=['start','stop','direction','taxon_id','source_y','version', 'entry_id', 'partial'])
        find_HDK(Cog_GK, args.source, args.nifgenes, args.gene, args.taxalevel, args.outfile, args.printcontigs, False, False)
    if args.source == "KeggGhostKoala":
        run_gene(Cog_GK, args.gene, args.source, args.specifictaxa, args.taxalevel, args.printcontigs, args.outfile)
    if args.source == "COG":
        print(CYEL + "\nNote that nifK is the beta chain of the nitrogenase and nifD is the alpha. Note COG2710 contains them both as 'Nitrogenase molybdenum-iron protein, alpha and beta chains'" + CEND)
        run_gene(Cog_GK, args.gene, args.source, args.specifictaxa, args.taxalevel, args.printcontigs, args.outfile)
    clean_up("Cog_GK.txt")

def clean_up(outfile):
    os.remove(outfile)
def check_files_existance(indirect):
    files = os.listdir(indirect)
    if not all(("genes_in_contigs.txt","genes_taxonomy.txt","gene_functions.txt","hmm_hits.txt","taxon_names.txt")) in files:
        pass
    else:
        print(CRED + "\nYou are missing a file that should have been exported from the Anvi'o database!\n")
        print("Make sure you have all of these files: genes_in_contigs.txt, genes_taxonomy.txt, gene_functions.txt, hmm_hits.txt, taxon_names.txt\n" + CEND)
        exit()

def check_source_input(data_frame, source):
    if not any(re.findall(r'|'.join(data_frame['source_x'].unique().tolist()+['KeggGhostKoala','COG']), source)):
        print(CRED + "\nPick from valid sources:")
        print(data_frame['source_x'].unique().tolist()+['KeggGhostKoala','COG'])
        exit()
    else:
        pass
    if source == "list":
        print('\033[93m'+ "\nPick from valid sources (case sensitive):")
        print((new3['source_x'].unique().tolist()+['KeggGhostKoala','COG'])) #lists sources to pick from
        print('\033[0m')
def check_foamlevel_input(Foamlevel):
    if Foamlevel == "list":
        print(CYEL + "\nLevel 2 of Foam assignments found in contigs where level 1 is Nitrogen cycle:")
        print(df_foam_hmms3['L2'].unique().tolist())
        exit()
    if Foamlevel is not None:
        df_foam_hmms4['L2'] = df_foam_hmms4.L2.str.replace("[(),]", " ")
        if not any(re.findall(r'|'.join(df_foam_hmms4['L2'].unique().tolist()), args.Foamlevel, re.IGNORECASE)):
            print(CRED + "\nPick from valid levels:")
            print(df_foam_hmms3['L2'].unique().tolist())
            exit()
    else:
        pass
def check_taxa_input(taxalevel):
    if taxalevel is not None:
        if any(re.findall(r'kingdom|phylum|class|order|family|genus|species', taxalevel, re.IGNORECASE)):
            pass
        else:
            print(CRED + "\nWhoops this doesn't look like a valid taxonomic level. Try one of these:")
            print('kingdom, phylum, class, order, family, genus, species'+ CEND)
            exit()
    else:
        pass
def check_gene_input(data_frame, source, gene):
    if source == "Nif_Hmms":
        if gene is None:
            pass
        elif gene is not None and not any(re.findall(r'|'.join(data_frame['Gene_description'].astype(str).str[0:4]),gene, re.IGNORECASE)) or gene == "list":
            genes = data_frame['Gene_description'].astype(str).str[0:4].unique().tolist()
            genes = [w.replace('tigr', 'recA') for w in genes]
            if gene == "list":
                print(CRED +"\nHere is list to pick from:")
            else:
                print(CRED +"\nGene not recognized :( \nHere is list to pick from:\n" + CEND)
            print(genes) #lists genes to pick from
            exit()
    if source == "Foam_Nitro":
        if gene is None:
            pass
        elif gene is not None and not any(re.findall(r'|'.join(data_frame['Gene'].unique().tolist()), gene, re.IGNORECASE)) or gene == "list":
            if gene == "list":
                print(CRED +"\nHere is list to pick from:")
            else:
                print(CRED +"\nGene not recognized :( \nHere is list to pick from:")
            print(data_frame['Gene'].unique().tolist())
            exit()

def print_contigs(data_frame, printcontigs, outfile):
    if printcontigs == True:
        if isinstance(data_frame, pd.DataFrame): #check that dataframe is passed
            if 'contig' in data_frame.columns:    #check to see if dataframe has contig column
                if outfile is None: #check that outfile name was giv
                    print(CRED +"\nYou need to select and out file name to write contigs to!\n"+ CEND)
                    exit()
                if outfile is not None:
                    with open('contigs_' + outfile, "w") as text_file: #writing out list of contigs that contain nif genes
                        pd.DataFrame(data_frame.contig.unique().tolist()).to_csv('contigs_'+outfile, header=0, sep='\t', mode='w', index=None)
                        print(CYEL + "\nCool we just printed some contigs into the file contigs_"+outfile +" . Now you can search for these contigs in MAG bins!")
            else:
                print(CRED +"There is no contig column in dataframe"+ CEND)
        else:
            print(CRED +"Yikes, this isn't a dataframe :("+ CEND)
    else:
        pass
def print_data_frame(data_frame, outfile):
    print(data_frame.head(10))
    Q = input(CPUR +'\nIs this the data frame you want printed to a csv file (yes/no)?' + CEND).lower()
    if Q.startswith('y'):#check to see if this is the dataframe you want printed
        if outfile is None: #check that outfile name was giv
            print(CRED +"\n You need to select an out file name to write dataframe to!\n"+ CEND)
            exit()
        if outfile is not None:
            with open(outfile, "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(data_frame.to_csv(outfile, header=0, sep='\t', mode='w', index=None))
                print(CYEL + "\nCool we just printed a data frame of contigs with their assigned functions and taxonmy into the file "+ outfile +CEND)
    if Q.startswith('n'):
        Q2 = input(CPUR +'\nDo you want to see the whole thing (yes/no)?' + CEND).lower()
        if Q2.startswith('y'):
            print(data_frame)
        if Q2.startswith('n'):
            print('\033[95m'+ "Sorry about that let's try again." + CEND)
def print_exit(data_frame, printcontigs, outfile, taxalevel):
    print_data_frame(data_frame, outfile)
    print_contigs(data_frame, printcontigs, outfile)
    if taxalevel is None:
        exit()
    else:
        pass
def print_taxa(data_frame, taxalevel, outfile):
    Q = input(CPUR +'\nDo you want unique taxa with hits to choosen gene(s) (option 1) or counts of how many times the taxa is assigned to a contig (option 2)?' + CEND).lower()
    if Q.startswith('2'):
        print('\033[93m' + "\nThis is a list of unique taxa with hits in your choosen gene(s) at this taxonomic level:\n" + CEND)
        print(data_frame['t_'+taxalevel].value_counts())
    if Q.startswith('1'):
        print('\033[93m' + "\nThis is a list of unique taxa with hits in your choosen gene(s) at this taxonomic level:\n" + CEND)
        print(data_frame['t_'+taxalevel].unique().tolist()) #get list of unique taxa that contain nifH.
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
def find_HDK(data_frame, source, nifgenes, gene, taxalevel, outfile, printcontigs, quite, getcontigs):
    if nifgenes == True: # This will find contigs with the same taxa that have hits in HDK
        if gene is not None:
            print(CRED + "\n You already said you wanted to look for nifHDK so i'm ignoring your request for another gene" + CEND)
        if taxalevel is not None:
            if source == "Nif_Hmms":
                nifH =  data_frame.loc[data_frame['Gene_description'].str.contains('NifH', flags=re.IGNORECASE) == True]
                nifK =  data_frame.loc[data_frame['Gene_description'].str.contains('nifK', flags=re.IGNORECASE) == True]
                nifD =  data_frame.loc[data_frame['Gene_description'].str.contains('nifD', flags=re.IGNORECASE) == True]
                nifHD = pd.merge(nifH, nifD, on=['contig'], how='inner')
                nifHDK = pd.merge(nifHD, nifK, on=['contig'], how='inner')
            if source == "Foam_Nitro":
                nifH =  data_frame.loc[data_frame['Gene'].str.contains('nifH', flags=re.IGNORECASE) == True]
                nifK =  data_frame.loc[data_frame['Gene'].str.contains('nifK', flags=re.IGNORECASE) == True]
                nifD =  data_frame.loc[data_frame['Gene'].str.contains('nifD', flags=re.IGNORECASE) == True]
                nifHD = pd.merge(nifH, nifD, on=['contig'], how='inner')
                nifHDK = pd.merge(nifHD, nifK, on=['contig'], how='inner')
            if source == "KeggGhostKoala":
                nifH =  data_frame.loc[data_frame['function'].str.contains('nifH', flags=re.IGNORECASE) == True]
                nifK =  data_frame.loc[data_frame['function'].str.contains('nifK', flags=re.IGNORECASE) == True]
                nifD =  data_frame.loc[data_frame['function'].str.contains('nifD', flags=re.IGNORECASE) == True]
                nifHD = pd.merge(nifH, nifD, on=['contig'], how='inner')
                nifHDK = pd.merge(nifHD, nifK, on=['contig'], how='inner')
            if source == "COG":
                nifH =  data_frame.loc[data_frame['function'].str.contains('nifH', flags=re.IGNORECASE) == True]
                nifDK =  data_frame.loc[data_frame['accession'].str.contains('COG2710', flags=re.IGNORECASE) == True]
                nifHDK = pd.merge(nifH, nifDK, on=['contig'], how='inner')
            if quite == True:
                taxa_list = nifHDK['t_'+taxalevel+'_x'].unique().tolist()
                return taxa_list
            if getcontigs == True:
                contigs_list = nifHDK['contig'].unique().tolist()
                return contigs_list
            if quite == False and getcontigs == False:
                print('\033[93m'+ "\nThe unique taxa at the level of " + taxalevel + " that contain nifHDK based on " + source + " assignments are: ")
                print(nifHDK['t_'+taxalevel+'_x'].unique().tolist())
                print('\n' + CEND)
                Q = input(CPUR +'\nDo you want to see the data (yes/no)?' + CEND).lower()
                if Q.startswith('y'):#check to see if this is the dataframe you want printed
                    print_exit(nifHDK, printcontigs, outfile, taxalevel)
                    exit()
                if Q.startswith('n'):
                    exit()
                if outfile is None: #check that outfile name was given
                    print(CRED +"\n You need to select an out file name to write dataframe to!\n"+ CEND)
                    exit()
                if outfile is not None:
                    with open(outfile, "w") as text_file: #writing out list of contigs that contain nif genes
                        pd.DataFrame(data_frame.to_csv(outfile, header=True, sep='\t', mode='w', index=None))
                    print(CYEL + "\nCool we just printed your data frame to the file "+outfile +CEND)
                    exit()
        if taxalevel is None:
            print(CRED + "\nYou need to pass a taxa with -t when you use --find-nifHDK-in-contigs")
            exit()
    else:
        pass

def specific_taxa(data_frame, specifictaxa, taxalevel):
    if specifictaxa is None:
        pass
    elif taxalevel is None:
        print(CRED +"\nYou need to pass a taxalevel when requesting a specific taxa")
        exit()
    elif specifictaxa is not None:
        print(CYEL + "\nThis is a dataframe, of contigs that had hits with from choosen gene(s) hmms hits for " + specifictaxa + "./n"+ CEND)
        data_frame = data_frame[data_frame['t_'+taxalevel].str.contains(specifictaxa)]
        print_data_frame(data_frame)
        exit()
def run_gene(data_frame, gene, source, specifictaxa, taxalevel, printcontigs, outfile):
    if gene is not None:
        if source == "KeggGhostKoala" or source == "COG":
            df = data_frame[data_frame['function'].str.contains(gene, flags=re.IGNORECASE)]
        if source == "Nif_Hmms":
            df = data_frame[data_frame['Gene_description'].str.contains(gene, flags=re.IGNORECASE)]
        if source == "Foam_Nitro":
            df = data_frame[data_frame['Gene'].str.contains(gene, flags=re.IGNORECASE)] #searching for particular rows with matching gene_name
        specific_taxa(df, specifictaxa, taxalevel)
        print(CYEL + "\nThis is a dataframe, of contigs that had hits with from " + source + " for your choosen gene(s). \nTo find out what taxa the contigs belong to pass -t with a taxalevel. If you already passed -t just hang on and we will get there\n" + CEND)
        print_exit(df, printcontigs, outfile, taxalevel)
        print_taxa(df, taxalevel, outfile)
    if gene is None:
        specific_taxa(data_frame, specifictaxa, taxalevel)
        print(CYEL + "\nThis will return a full dataframe, of contigs that had hits with from " + source + ".\nTo reduce the size of the dataframe pick a gene and/or taxalevel to look at.\n If you already passed -t just hang on and we will get there.\n" + CEND)
        print_exit(data_frame, printcontigs, outfile, taxalevel)
        print_taxa(data_frame, taxalevel, outfile)

if __name__ == '__main__':
    main()
