#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019
## Given a Kaiju .summary file script will give a .csv file to read into R

#importing packages
import glob,os
import pandas as pd
import re
import sys
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="exsum.py",description="Parses Kaiju .summary files into an R readable .csv file")
    parser.add_argument("-t",'--taxonomic-level', dest="taxlevel", action="store", required=True, help='Pick One: phylum, class, order, family, genus, species (required)')
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were .summary files are found (required)")
    return parser.parse_args()

def input_output_dir(indirname):
    #Have we got an input and output directory? If not, exit.
    if indirname is None:
        print("No input directory name (exiting)")
        sys.exit(1)
def call_in_files(taxlevel,indirname):
    #getting name of file from taxlevel
    line = '*_kaiju_out_.summary'
    index = line.find('.summary')
    output_line = line[:index] + taxlevel + line[index:]
    #Calling in files
    os.chdir(indirname)
    global infiles
    infiles = glob.glob(output_line)

def main():
    #Parse command-line
    args = parse_cmdline()
    input_output_dir(args.indirname) #Have we got an input and output directory? If not, exit.
    call_in_files(args.taxlevel, args.indirname) #call in files
    summary = [] #making blank dataframe to store counts of each gene from output of for loop
    for file in infiles:
        #filename = f.name + '.csv'
        df = pd.read_csv(file, sep='\t', header=0)
        df.columns = ["Sample","Percent_Reads", "Reads", "Taxon_id", args.taxlevel.capitalize()]
        df['Sample'] = df['Sample'].str.replace('_S(.*?)_kaiju_out.tsv', '') #editing "Sample_Names"
        #df['Taxanomic_Level'] = f.name
        df['Taxanomic_Level'] = args.taxlevel.capitalize()
        #Numbers at beginning of file are Month-Farm-Sample
        df['Farm'] = pd.np.where(df.Sample.str.startswith("8"), "Farm_8", #Adding column that will have the farm name
                           pd.np.where(df.Sample.str.startswith("6"), "Farm_6",
                           pd.np.where(df.Sample.str.startswith("5"), "Farm_5",
                           pd.np.where(df.Sample.str.startswith("1"), "Farm_1", "No_name"))))
        summary.append(df) #store DataFrame in list
        name = args.taxlevel + '_summary.csv'
        print("\n"+'\033[93m' + "Finished parsing " + open(file,'r').name + "\n")
    summary = pd.concat(summary)
    #write DataFrame to comma separated file (.csv) with file name and TIGRFAM counts
    summary.to_csv(name, sep=',' ,index=False)
    print("\n"+'\033[93m' + "New file created: " + name + "\n")

if __name__ == '__main__':
    main()
