#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019
## Given a Prokka .cog file script will give a .feather file to read into R

## .cog file created by the following line:
#egrep "COG[0-9]{4}" ../10_6_3_Output/10_6_3.gff | cut -f9 | sed 's/.\+COG\([0-9]\+\).\+;locus_tag=\(GIGOMEAJ_[0-9]\+\);.\+/\2\tCOG\1/g' > 10_6_3.cog
#You will need to replace 10_6_3 with your file name and "GIGOMEAJ" with whatever string prokka made a the beginning of your gene calls
#You can also make .cog file anyway you want. Just needs to have these columns:
#Contig_name \t cog####

#importing packages
import feather
import glob,os
import pandas as pd
from pandas import DataFrame
import re
import csv
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="Prokka_cog_counts.py",description="Parses Prokka .cog files into an R readable .feather file")
    parser.add_argument("-c", dest="cats", action="store", required=True, help='Give full path to Cog_Cats.txt file. This is a tab separated file (required)')
    parser.add_argument("-cg", "--cognames", dest="cognames", action="store", default=None, required=True, help="Give full path to cognames2003-2014.tab file (required)")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were .cog files are found (required)")
    return parser.parse_args()

def split_it(x):
    return re.sub(r'([A-Z])([A-Z])', r'\1,\2',x)
    return re.sub(r'([A-Z])([A-Z])([A-Z])', r'\1,\2,\3',x)

if __name__ == '__main__':
    #Parse command-line
    args = parse_cmdline()
    #making dataframes from csv files
    df_cat = pd.read_csv(args.cats, sep='\t', header=0)
    df_let = pd.read_csv(args.cognames, sep='\t', header=0, encoding='windows-1252')
    #We need to massage the df_let file a bit changing column names
    df_let.columns = ['COG_ID','Letter','Name'] #rename rows
    #Adding addtional rows for cases where there are two letters in letter column
    df_let['Letter'] = df_let['Letter'].apply(split_it) ##first we separate multiple capitial letters in the "Letters" columns by commas by applying a function across the dataframe
    new_let = DataFrame(df_let.Letter.str.split(',').tolist(), index=([df_let.Name,df_let.COG_ID])).stack()
    new_let = new_let.reset_index()[[0, 'COG_ID', 'Name']] # "Letter" variable is currently labeled 0
    new_let.columns = ['Letter','COG ID','Name'] # renaming Letter
    #Now we can move one to combining dataframe
    df = pd.merge(df_cat, new_let, on='Letter', how="right")
    df.to_csv('cogs_and_cats.csv', sep=',' ,index=False) # writing to a .csv file
    #Calling in files
    os.chdir(args.indirname) #change directory
    infiles = glob.glob("*.cog")
    for file in infiles:
        df_cogs = pd.read_csv(file, sep='\t', header=None)
        #changing column names
        df_cogs.columns = ['Gene_Call','COG ID']
        df_cogs_new = pd.merge(df, df_cogs, on='COG ID', how="right")
        df_cogs_new.to_csv('genecalls_to_cats.csv', sep=',' ,index=False) # writing to a .csv file
        #Note that df_cogs_new will have a different shape than df_cogs since some cogs have multiple categories assigned to them
        #getting frequency counts
        counts = pd.DataFrame(df['Catergory'].value_counts())
        counts.reset_index(inplace=True)
        counts.columns = ['Catergory','Counts']
        counts['Percent'] =  (counts['Counts']/sum(counts['Counts']))*100
        #Write to file for reading from R
        feather.write_dataframe(counts, "Prokka_cog_counts.feather")
