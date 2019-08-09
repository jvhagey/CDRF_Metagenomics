#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/

#In shell terminal grep to get name and accession numbers of hmms
#grep -A1 "NAME" FOAM-hmm_rel1a.hmm > Acc_num.txt

#importing packages
import glob,os
import pandas as pd
import re
from argparse import ArgumentParser


#Defining function to process command-line arguments
def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="FOAM_to_anvio.py", description="Converts Acc_num file into genes.txt file that Anvi'o can use")
	parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were .fasta files are found (required)")
	return parser.parse_args()

if __name__ == '__main__':
	#Parse command-line
	args = parse_cmdline()
	#Calling in files
	os.chdir(args.indirname) #change directory
	#making txt file for gene.txt file to run anvi-run-hmm
	file = "Acc_num.txt"
	with open(file, "r+") as f: #open file for reading and writing
		filedata = f.read()
		#Removing lines starting with --
		filedata = re.sub(r'--\n', '', filedata, flags=re.MULTILINE)
		#replacing two or three spaces with only one
		filedata = re.sub(r'   ', ' ', filedata, flags=re.MULTILINE)
		filedata = re.sub(r'  ', ' ', filedata, flags=re.MULTILINE)
	with open(file, 'w') as f2: #writing new header for edited file removing blank lines at beginning and end
		f2.writelines(filedata)
	df = pd.read_csv(file, sep=' ', header=None)
	df.columns = ['Columns','Rows'] #rename columns
	df = df.pivot(columns = "Columns", values = "Rows")
	#remove none values and move up the cells
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')
	df = df[['NAME', 'ACC']]
	df["hmmsource"] = "FOAM_DB"
	df.columns = ["gene","accession","hmmsource"] #rename columns
	#Anvi'o won't like duplicate names in the "gene" column so we will have to add the accession to the end of each KO
	df['gene'] = df.gene.map(str) + "_" + df.accession
	#write DataFrame to tab separated file (.csv)
	df.to_csv('gene_FOAM.txt', sep='\t',index=False)
