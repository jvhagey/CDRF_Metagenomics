#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019
## given an hmm file with repeated NAMES this script outputs a new hmm file with NAMES numbered

#importing packages
import os
import re
import csv
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="Identify_Nif_Taxonomy.py", description="Parses Hmmscan .tblout files into an R readable .feather file")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were at .txt files are found (required)")
    parser.add_argument("-f", "--infile", dest="infile", action="store", default=None, required=True, help="Input file name were at .hmm file is found (required)")
    parser.add_argument("-o", "--outfile", dest="outfile", action="store", default=None, required=True, help="out file name for new .hmm file with names changed (required)")
    return parser.parse_args()

if __name__ == '__main__':
    #Parse command-line
    args = parse_cmdline()
    #Calling in files
    os.chdir(args.indirname) #change directory
    file = args.infile
    open(args.outfile, 'a').close() #making empty text file for later
    with open(file) as f:
        seen = set()
        counting = list()
        for line in f:
            if line.startswith('NAME'):# grab name line
                line = re.sub(r'\n', '', line, flags=re.MULTILINE) #grab name line remove \n at end of each line
                if line in seen:
                    counting.append(line) #add line to counting list
                    num = counting.count(line) #count the number of times you see line in counting list
                    new_line = line + "_" + str(num) #make new name with the number on end
                    with open("genes_test_new.hmm", "a") as f2:
                #write line to output file
                        f2.write(new_line)
                        f2.write("\n")
                else:
                    seen.add(line)
                    counting.append(line)
                    new_line1 = line + "_1" #make new name with the number on end
                    with open("genes_test_new.hmm", "a") as f3:
                        f3.write(new_line1)
                        f3.write("\n")
            else:
                with open("genes_test_new.hmm", "a") as f4:
                    f4.write(line)
