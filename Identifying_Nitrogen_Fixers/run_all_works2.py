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
from contextlib import contextmanager

def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="Find_common_contigs.py", description="Parses anvi-export-table table files")
    parser.add_argument("-i", "--indir", dest="indirname",action="store", default=None, required=True, help="Input directory name were taxon_names, genes_taxonomy and hmm_hits .txt files are found (required)")
    parser.add_argument("-o", "--outfile", dest="outfile",action="store", default=None, required=False, help="File name for output to be written to.")
    parser.add_argument("-t", "--taxa-level",type=str, dest="taxalevel", action="store", default=None, help='Pick taxa level you want reported')
    parser.add_argument("--get-HDK-contigs", dest='getHDKcontig',action="store_true", default=False, help="Get contigs with nifHDK genes all IN THE SAME CONTIG and print to text file.")
    parser.add_argument("--get-HDK-taxa", dest='getHDKtaxa',action="store_true", default=False, help="Same as the --get-HDK-contigs, but you get the option of seeing taxa with nifHDK genes IN THE SAME CONTIG as well.")
    parser.add_argument("--print-HDK-data-frame", dest='printdataframe',action="store_true", default=False, help="Print dataframe with contigs that contain nifHDK genes IN THE SAME CONTIG at chosen taxa level \
                                                                                                                    to text file. This will print the dataframe that --get-HDK-contigs is based on.")
    parser.add_argument("--print-common-taxa", dest='printtaxa',action="store_true", default=False, help="Find contigs with nifHDK genes IN SAME CONTIG, identify taxa assiged to these gene calls at level chosen.\
                                                                                                            Then find common to find the taxa that are identified with all 4 annotation methods. Does take the output of --get-HDK-taxa.")
    parser.add_argument("--print-taxa-nifHDK-diff-contigs", dest='taxamanygenes',action="store_true", default=False, help="Get taxa for nifHDK genes that are found in different contigs (Honestly, you can only do it for 2 genes).")
    parser.add_argument("-g", "--gene", type=str, dest="gene", action="store", default=None, required=False, help="Name of gene to search in dataframe for. You can pass multiple genes by separating them by a pipe. Ex: 'geneA|geneB'")
    parser.add_argument("--print-gene-contigs", dest='printgenecontigs',action="store_true", default=False, help="Get contigs containing gene(s) of choice to text file.")
    args = parser.parse_args()
    return args

#set colors for warnings so they are noticable
CRED = '\033[91m'+'\nWarning:'
CYEL = '\033[93m'
CPUR = '\033[95m'+'\nQuick Question:'
CEND = '\033[0m'

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
def shorten(s, subs):
    i = s.index(subs)
    return s[:i+len(subs)]

def arg_check(printtaxa, getHDKcontig, gene, taxalevel, printdataframe, getHDKtaxa, taxamanygenes):
    if gene == "nifK" or gene == "nifD":
        print(CYEL+ '\n Just a heads up that nifK and nifD are not found in the COGs. COG2710 You need to pass "nifK|COG2710" to get what you want! \n'+ CEND)
        exit()
    if getHDKtaxa == True and taxalevel is None:
        print(CRED+ "\n You need to pick a taxa level by passing -t when choosing --get-HDK-taxa, cuz ya know I won't know what to print otherwise. ¯\_(ツ)_/¯ \n"+ CEND)
        exit()
    if printtaxa == True and taxalevel is None:
        print(CRED+ '\n You need to pick a taxa level by passing -t when choosing --print-taxa mode \n'+ CEND)
        exit()
    if printtaxa == True and getHDKcontig == True:
        print(CRED + "\n You can't pass both --get-taxa-contigs and --print-taxa. Choose just one at a time :( \n"+ CEND)
        exit()
    if getHDKtaxa == True and taxalevel == False:
        print(CRED + "\n You need pass both --get-HDK-taxa and a taxa level with -t.\n"+ CEND)
        exit()
    if printtaxa == False and getHDKcontig == False and gene is None and printdataframe == False and getHDKtaxa == False and taxamanygenes == False:
        print(CRED + "\n You need to pick a mode, either pass gene name with '-g gene_x' or pick one of the following modes --get-HDK-taxa, --get-HDK-contigs, --print-gene-contigs , --print-HDK-data-frame, or --get-contigs or --print-taxa. Choose just one at a time :( \n"+ CEND)
        exit()
    if taxamanygenes == True:
        if taxalevel is None:
            print(CRED + "\n You need to pick a taxa level and pass with -t. \n"+ CEND)
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
        global df_nif_hmms2, df_foam_hmms4, df_foam_hmms2, GK, Cog
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
def get_gene_df(gene, taxalevel, printgenecontigs):
    global Foam_gene_list, GK_gene_list, Cog_gene_list, Nif_Hmms_gene_list, Foam_gene, GK_gene, Cog_gene, Nif_Hmms_gene
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
        if "COG" in gene:
            Cog_gene =  Cog.loc[Cog['accession'].str.contains(gene, flags=re.IGNORECASE) == True]
            Cog_gene_list = Cog_gene.contig.unique().tolist()
        #print(CYEL + "\n This is a list of contigs found with Cogs for "+gene+":\n")
        #print(Cog_gene.contig.unique().tolist())
        #global Nif_Hmms_gene_list
        Nif_Hmms_gene = df_nif_hmms2[df_nif_hmms2['Gene_description'].str.contains(gene, flags=re.IGNORECASE)]
        Nif_Hmms_gene_list = Nif_Hmms_gene.contig.unique().tolist()
        #print(CYEL + "\n This is a list of contigs found with nif TigrFams for "+gene+":\n")
        #print(Nif_Hmms_gene.contig.unique().tolist())
        if taxalevel is not None:
            global Foam_taxa_list, GK_taxa_list, Cog_taxa_list, Nif_Hmms_taxa_list
            Nif_Hmms_taxa_list = Nif_Hmms_gene['t_'+taxalevel].unique().tolist()
            Cog_taxa_list = Cog_gene['t_'+taxalevel].unique().tolist()
            GK_taxa_list = GK_gene['t_'+taxalevel].unique().tolist()
            Foam_taxa_list = Foam_gene['t_'+taxalevel].unique().tolist()
        if printgenecontigs == True:
            with open('contigs_Nif_Hmms_'+gene+'.txt', "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(Nif_Hmms_gene.contig.unique().tolist().to_csv('contigs_Nif_Hmms_'+gene+'.txt', header=True, sep='\t', mode='w', index=None))
            with open('contigs_Cog_'+gene+'.txt', "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(Cog_gene.contig.unique().tolist().to_csv('contigs_Cog_'+gene+'.txt', header=True, sep='\t', mode='w', index=None))
            with open('contigs_GK_'+gene+'.txt', "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(GK_gene.contig.unique().tolist().to_csv('contigs_GK_'+gene+'.txt', header=True, sep='\t', mode='w', index=None))
            with open('contigs_Foam_Hmms_'+gene+'.txt', "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(Foam_gene.contig.unique().tolist().to_csv('contigs_Foam_Hmms_'+gene+'.txt', header=True, sep='\t', mode='w', index=None))
            print(CYEL + "\nCool we just printed some contigs into individual files. Now you can search for these contigs in MAG bins!")
            exit()
        else:
            pass
def get_taxa_or_contigs(taxalevel, printtaxa, getHDKcontig, printdataframe, getHDKtaxa):
    global Foam_list, GK_list, Cog_list, Nif_Hmms_list
    if printdataframe == True:
        Foam_list = Anvi_table_parser_new3.find_HDK(df_foam_hmms4, "Foam_Nitro", True, None, None, 'nifHDK_Foam_Nitro.txt', False, False, getHDKcontig, printdataframe, False)
        #global GK_list
        GK_list = Anvi_table_parser_new3.find_HDK(GK, "KeggGhostKoala", True, None, None, 'nifHDK_GK.txt', False, False, getHDKcontig, printdataframe, False)
        #global Cog_list
        Cog_list = Anvi_table_parser_new3.find_HDK(Cog, "COG", True, None, None, 'nifHDK_Cog.txt', False, False, getHDKcontig, printdataframe, False)
        #global Nif_Hmms_list
        Nif_Hmms_list = Anvi_table_parser_new3.find_HDK(df_nif_hmms2, "Nif_Hmms", True, None, None, 'nifHDK_nif_hmms.txt', False, False, getHDKcontig, printdataframe, False)
    else:
        Foam_list = Anvi_table_parser_new3.find_HDK(df_foam_hmms4, "Foam_Nitro", True, None, taxalevel, None, False, printtaxa, getHDKcontig, printdataframe, getHDKtaxa)
        #global GK_list
        GK_list = Anvi_table_parser_new3.find_HDK(GK, "KeggGhostKoala", True, None, taxalevel, None, False, printtaxa, getHDKcontig, printdataframe, getHDKtaxa)
        #global Cog_list
        Cog_list = Anvi_table_parser_new3.find_HDK(Cog, "COG", True, None, taxalevel, None, False, printtaxa, getHDKcontig, printdataframe, getHDKtaxa)
        #global Nif_Hmms_list
        Nif_Hmms_list = Anvi_table_parser_new3.find_HDK(df_nif_hmms2, "Nif_Hmms", True, None, taxalevel, None, False, printtaxa, getHDKcontig, printdataframe, getHDKtaxa)
def find_common(nif_hmms, foam_hmms, GK, cog, printtaxa, getHDKcontig, gene, outfile, taxalevel, quite):
    if printtaxa == True:
        word = "taxa"
    if getHDKcontig == True:
        word = "contigs"
    if gene is not None:
        word = "contigs for " + gene
    if gene is not None and taxalevel is not None:
        word = "taxa with " + gene
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
    global final_list
    final_list = set1.intersection(set2)
    print("\nThe " +word+ " that were found in common will all three methods are:")
    print(list(final_list))
    if quite == False:
        Q = input(CPUR +'\nDo you want to print the contigs (yes/no)?' + CEND).lower()
        if Q.startswith('y'):#check to see if this is the dataframe you want printed
            if outfile is not None:
                final_df = pd.DataFrame(list(final_list), index=None)
                with open(outfile, "w") as text_file: #writing out list of contigs that contain nif genes
                    pd.DataFrame(final_df.to_csv(outfile, header=True, sep='\t', mode='w', index=None))
                    print(CYEL + "\nCool we just printed your contigs to the file "+outfile +CEND)
            if outfile is None:
                print(CRED +'\n You need to pass a outfile name with the -o option \n' + CEND)
                exit()
            pass
        if Q.startswith('n'):
            pass
def find_HDK_taxa(taxalevel, outfile):
    df = df_nif_hmms2[df_nif_hmms2['contig'].str.match(r'|'.join(final_list), flags=re.IGNORECASE)]
    df = df.drop(columns=['gene_callers_id','taxon_id','t_class','t_order','entry_id','source_x','gene_hmm_id','gene_name'])
    print('\n')
    print(CYEL+'\n This is a dataframe of contigs that have hits with HDK within the same contig')
    print(df)
    print(CYEL+'\n These are the unique ' + taxalevel + ' that have HDK genes within the same contig')
    print(df['t_'+taxalevel].unique().tolist())
    Q = input(CPUR +'\nDo you want to print the dataframe (yes/no)?' + CEND).lower()
    if Q.startswith('y'):#check to see if this is the dataframe you want printed
        if outfile is not None:
            with open(outfile, "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(df.to_csv(outfile, header=True, sep='\t', mode='w', index=None))
                print(CYEL + "\nCool we just printed your dataframe to the file "+outfile +CEND)
        if outfile is None:
            print(CRED + '\n You need to pass a outfile name with the -o option \n'+ CEND)
            exit()
        pass
    if Q.startswith('n'):
        pass
def run_gene_workflow(gene, taxalevel, outfile):
    if gene is not None:
        Q = input(CPUR +'\nDo you want to contigs containing your gene or taxa with this gene (taxa/contigs)?' + CEND).lower()
        if Q.startswith('c'):#check to see if this is the dataframe you want printed
            find_common(Nif_Hmms_gene_list, Foam_gene_list, GK_gene_list, Cog_gene_list, False, False, gene, outfile, None)
        if Q.startswith('t'):
            find_common(Nif_Hmms_taxa_list, Foam_taxa_list, GK_taxa_list, Cog_taxa_list, False, False, gene, outfile, taxalevel, False)
            Q = input(CPUR +'\nDo you want to counts of taxa containing your gene (yes/no)?' + CEND).lower()
            if Q.startswith('y'):
                print(CYEL+"\nThe " +taxalevel+ " that have " +gene+ " based on custom nif Hmms:\n")
                print(Nif_Hmms_gene['t_'+taxalevel].value_counts().nlargest(8))
                print(CYEL+"\nThe " +taxalevel+ " that have " +gene+ " based on Cogs:\n")
                print(Cog_gene['t_'+taxalevel].value_counts().nlargest(8))
                print(CYEL+"\nThe " +taxalevel+ " that have " +gene+ " based on GhoatKOALA:\n")
                print(GK_gene['t_'+taxalevel].value_counts().nlargest(8))
                print(CYEL+"\nThe " +taxalevel+ " that have " +gene+ " based on Foam Hmms:\n")
                print(Foam_gene['t_'+taxalevel].value_counts().nlargest(8))
                print(CYEL+"\nThe " +taxalevel+ " that have " +gene+ " common to all 4 methods:\n")
                save = pd.DataFrame(Foam_gene['t_'+taxalevel].value_counts())
                save.index.name = taxalevel.capitalize()
                save.reset_index(inplace=True)
                save.columns = [taxalevel.capitalize(),"Counts"]
                with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # forcing pandas to print all output
                    print(pd.DataFrame(save.loc[save[taxalevel.capitalize()].str.contains("|".join(final_list))]))
            if Q.startswith('n'):
                exit()
def find_HDK_taxa_diff_contigs(taxalevel, outfile):
    get_gene_df("nifH", taxalevel, False)
    with suppress_stdout():
        find_common(Nif_Hmms_taxa_list, Foam_taxa_list, GK_taxa_list, Cog_taxa_list, False, False, "nifH", outfile, taxalevel, True)
    save = pd.DataFrame(Foam_gene['t_'+taxalevel].value_counts())
    save.index.name = taxalevel.capitalize()
    save.reset_index(inplace=True)
    save.columns = [taxalevel.capitalize(),"Counts"]
    taxa_with_gene_nifH = pd.DataFrame(save.loc[save[taxalevel.capitalize()].str.contains("|".join(final_list))])
    get_gene_df("nifK|COG2710|nifD", taxalevel, False)
    with suppress_stdout():
        find_common(Nif_Hmms_taxa_list, Foam_taxa_list, GK_taxa_list, Cog_taxa_list, False, False, "nifK|COG2710|nifD", outfile, taxalevel, True)
    taxa_with_gene_nifKD = pd.DataFrame(save.loc[save[taxalevel.capitalize()].str.contains("|".join(final_list))])
    print(CYEL + "\nThese are the "+taxalevel+" that have contigs with nifH, nifD and nifK in multiple contigs. These were found in common will all 4 annotation methods. \n")
    print(set(taxa_with_gene_nifH[taxalevel.capitalize()]).intersection(set(taxa_with_gene_nifKD[taxalevel.capitalize()])))
    print("\n")
    exit()

def get_gene_call(dataframe, outfile):
    df = dataframe[dataframe['contig'].str.contains(r'|'.join(final_list))]
    #df['gene_callers_id'] = '>' + df['gene_callers_id'].astype(str) #adding ">" at the start of gene call so that it can be a searchable fasta file
    Q = input(CPUR +'\nDo you want the whole dataframe or one dataframe per gene (whole/per gene)?' + CEND).lower()
    if Q.startswith('p'):#check to see if this is the dataframe you want printed
        list = df["Gene_description"].unique().tolist()
        y = []
        for lst in list:
            z = lst.split(":",1 )[0]
            z = z.replace(":", "")
            y.append(z)
            testList = []
            testList.append(y)
        for i in y:
            df2 = df[df["Gene_description"].str.contains(i)]
            df2 = df2["gene_callers_id"]
            file_name = i +"_genecalls.txt"
            with open(file_name, "w") as text_file: #writing out list of contigs that contain nif genes
                pd.DataFrame(df2.to_csv(file_name, header=0, sep='\n', mode='w', index=None))
                print(CYEL + "\nPrinted a gene calls for " + i + " into the file "+ file_name +CEND)
        exit()
    if Q.startswith('w'):#check to see if this is the dataframe you want printed
        print(CYEL +"\nHere is part of the dataframe with contigs and gene calls.\n")
        print(CYEL + "\nNote that printing this will over ride file that you have printed contigs.\n")
        Anvi_table_parser_new3.print_data_frame(df, outfile)
        exit()

def main():
    args = parse_cmdline()
    os.chdir(args.indirname) #change directory
    arg_check(args.printtaxa, args.getHDKcontig, args.gene, args.taxalevel, args.printdataframe, args.getHDKtaxa, args.taxamanygenes)
    get_datatables(args.indirname)
    Anvi_table_parser_new3.check_gene_input(df_foam_hmms4, "Foam_Nitro", args.gene)
    Anvi_table_parser_new3.check_gene_input(df_nif_hmms2, "Nif_Hmms", args.gene)
    if args.taxamanygenes == True:# NOTE:
        find_HDK_taxa_diff_contigs(args.taxalevel, args.outfile)
    get_gene_df(args.gene, args.taxalevel, args.printgenecontigs)
    run_gene_workflow(args.gene, args.taxalevel, args.outfile)
    if args.printtaxa == True or args.getHDKcontig == True:
        get_taxa_or_contigs(args.taxalevel, args.printtaxa, args.getHDKcontig, args.printdataframe, False)
        find_common(Nif_Hmms_list, Foam_list, GK_list, Cog_list,args.printtaxa, args.getHDKcontig, None, None, None, False)
        get_gene_call(df_nif_hmms2, args.outfile)
        get_gene_call(df_foam_hmms4, args.outfile)
    if args.printdataframe == True:
        get_taxa_or_contigs(args.taxalevel, args.printtaxa, args.getHDKcontig, args.printdataframe, args.getHDKtaxa)
    if args.getHDKtaxa == True:
        get_taxa_or_contigs(False, False, True, False, args.getHDKtaxa)
        find_common(Nif_Hmms_list, Foam_list, GK_list, Cog_list, False, True, None, args.outfile, None, False)
        find_HDK_taxa(args.taxalevel, args.outfile)

if __name__ == '__main__':
    main()
