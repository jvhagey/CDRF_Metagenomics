#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2020
## Given a hmmr .tblout file from searching against FOAM HMMs, script will give a .feather file to read into R
#Script takes tblout output file from hmmscan and converts them to a table that displays gene counts that can be transfered into R for analysis/graphing

#importing packages
import glob,os
import pandas as pd
import re
import csv
#import feather
import sys
import numpy as np
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="Kaiju_to_outtab.py",description="Parses kaiju-names_out files into an R readable .feather file, makes files otu_table.txt and tax_table.txt to make phyloseq object. NOTE: FUNGI ARE CURRENTLY DELETED FROM LIST")
    parser.add_argument("-f", dest="filename", action="store", required=True, help='Give non-unique file name extension should be in the same folder.')
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were kaiju-names_out.tsv files are found (required)")
    #parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None, required=True, help="Output directory (required)")
    return parser.parse_args()
    return args

def clean_kaiju():
#Calling in files
  Kaiju_files = glob.glob("*_kaiju-names_out.tsv") #files with extention that will be looped through code
  #Kaiju_files = glob.glob("*" + args.filename) #files with extention that will be looped through code
  strip_character = "\t"
  temp = open("temp_file.txt", 'w')
  temp.close()
  for file in Kaiju_files:
    with open(file,"r") as infile, open("temp_file.txt", 'a') as outfile:
      for line in infile:
        if line.startswith("C"): #grab only classified lines
          line = re.sub(';\s', ';', line)#remove space between colon and taxa name
          go_away = strip_character.join(line.split(strip_character)[:2])#define waht I want to remove
          line = line.replace(go_away,"") #remove it
          line = re.sub('\t[0-9]', '', line)#cut tab at the beginnning
          line = re.sub(';', '\t', line)#cut tab at the beginnning
          outfile.write(line)

def tax_table():
  Kaiju_files = glob.glob("*_kaiju-names_out.tsv") #files with extention that will be looped through code
  #Kaiju_files = glob.glob("*"+args.filename) #files with extention that will be looped through code
  taxa_col = ("Organism\tKingdom\tGroup\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies")
  with open("tax_table.txt", 'w') as new:
    new.write(taxa_col)
    new.write("\n")
  Fungi = open("Fungi_taxa.txt", 'w')
  Fungi.close()
  for file in Kaiju_files:
    with open("temp_file.txt","r") as infile, open("tax_table.txt", 'a') as outfile, open("Fungi_taxa.txt", 'a') as outfileFungi:
      for line in infile:
        line = re.sub('^[0-9]+\t', '', line)#cut numbers at beginning
        line = re.sub('Proteobacteria\t',"Proteobacteria: ", line) #fix so protiobacteria stays in phylum column
        n = line.count("\t")
        n_start = line.count("\t")
        if n <= 8:
          line = line.strip()
        while n < 8:
          line += "\tNA"
          n = line.count("\t")
        if n >= 8:
          line = line.strip() + "\n"# get rid of trailing tabs
        if "Fungi" in line: #Fungi have a different number of columns so just setting them aside.
            outfileFungi.write(line)
        if not "Fungi" in line:
            outfile.write(line)
  df = pd.read_csv("tax_table.txt", header=0, sep='\t')
  df = df.drop_duplicates(keep='first')#drop duplicate rows
  df['RN'] = np.arange(len(df))#make new column with rumbers
  df["RN"] = 'Kai_' + df["RN"].astype(str) #add prefix
  df.index = df['RN'] #make into row names
  df.to_csv('tax_table_RN.txt', sep='\t',index=True) #write DataFrame to comma separated file (.csv) with file name and FOAM hmm counts
  df = df.drop('RN',1)
  df.to_csv('tax_table.txt', sep='\t',index=True, header=True, quotechar='"', quoting=csv.QUOTE_NONNUMERIC) #write DataFrame to comma separated file (.csv) with file name and FOAM hmm counts
  #df.write_dataframe(df, "tax_table.feather")

def otu_table():
  Kaiju_files = glob.glob("*_kaiju-names_out.tsv") #files with extention that will be looped through code
  #Kaiju_files = glob.glob("*"+args.filename) #files with extention that will be looped through code
  print(Kaiju_files)
  otu_col = ("Sample\tReads\tOrganism\tKingdom\tGroup\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies")
  with open("otu_table_first.txt", 'w') as new:
    new.write(otu_col)
    new.write("\n")
  for file in Kaiju_files:
    sam_name = os.path.splitext(file)[0]
    #sam_name = re.sub('_kaiju-names_out', '', sam_name)
    sam_name = re.sub("-names_out_kaiju", '', sam_name)
    with open("temp_file.txt","r") as infile, open("otu_table_first.txt", 'a') as outfile, open("Fungi_taxa.txt", 'a') as outfileFungi:
      for line in infile:
        line = sam_name + "\t" + line
        n = line.count("\t")
        n_start = line.count("\t")
        if n <= 9:
          line = line.strip()
        while n < 9:
          line += "\tNA"
          n = line.count("\t")
        if n >= 9:
          line = line.strip() + "\n"# get rid of trailing tabs
        if "Fungi" in line: #Fungi have a different number of columns so just setting them aside.
          outfileFungi.write(line)
        if not "Fungi" in line:
            outfile.write(line)
  df = pd.read_csv("otu_table_first.txt", header=0, sep='\t')
  df_taxa = pd.read_csv("tax_table_RN.txt", header=0, sep='\t')
  df_merge = pd.merge(df, df_taxa, how='inner', on=['Organism','Kingdom','Group','Phylum','Class','Order','Family','Genus','Species'])
  df_merge.to_csv('check.txt', sep='\t',index=True)
  df_merge = df_merge.drop_duplicates(keep='first')
  df_otu = df_merge[['Sample', 'Reads','RN']]
  df_otu_wide = df_otu.pivot(index='RN', columns='Sample', values='Reads')
  df_otu_wide.to_csv('otu_table.txt', sep='\t',index=True)
  #feather.write_dataframe(df_otu_wide, "otu_table.feather")

def main():
  #args = parse_cmdline()
  os.chdir("./")
  clean_kaiju()
  tax_table()
  otu_table()
  os.remove("temp_file.txt")

if __name__ == '__main__':
  main()
