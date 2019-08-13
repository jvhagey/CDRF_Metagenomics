#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/

#importing packages
#import feather
import pandas as pd
import re
import glob,os
import csv

from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="GK_Count_Frequency.py", description="GhostKOALA output files into a KO count to read into R")
	parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name GhostKOALA output files are found (required)")
	parser.add_argument("-c", "--count", dest="countoption", action="store", default=None, required=True, help="Select '-c taxa' to return taxa counts or '-c cat' or category counts (required)")
	return parser.parse_args()

if __name__ == '__main__':
	#Parse command-line
	args = parse_cmdline()
	#Calling in files
	os.chdir(args.indirname) #change directory
	#making dataframes
	df_ortho = pd.read_csv("KeggOrthology_Table1.txt", sep=',', header=0)
	df_userdf = pd.read_csv("user_ko_definition_all.txt", sep='\t', header=None)
	df_taxdf = pd.read_csv("user.out.top_all.txt", sep='\t', header=None)
	#df_kodf = pd.read_csv("user_ko_all_cleaned.txt", sep='\t', header=None)
	#adding column names
	df_userdf.columns = ['Query','accession','Definition','Score','Second best','Second Score'] #renaming columns
	df_userdf = df_userdf.drop(columns=['Definition','Score','Second best','Second Score']) #dropping columns
	df_taxdf.columns = ['Query','accession','Kingdom','Taxa','Lower_taxa','TaxaID','Score'] #renaming columns
	df_taxdf['Query'] = df_taxdf['Query'].map(lambda x: x.lstrip('user:')) #remove "user:" at beginning of each row in column "Query"
	df_taxdf = df_taxdf.drop(columns=['TaxaID','Score']) #droping columns
	new = pd.merge(df_userdf, df_ortho, on='accession')# This takes to much memory so have to preprocess
	# creating a empty bucket to save result
	del(df_userdf) #delete to save memory
	del(df_ortho)
	new2 = pd.merge(new, df_taxdf, on='Query')
	del(new) #delete to save memory
	del(df_taxdf) #delete to save memory
	if args.countoption is None:
		print("You need to pick either taxa or cat for counts")
	if args.countoption == 'cat':
		cat1 = pd.DataFrame(new['Category1'].value_counts())
		cat2 = pd.DataFrame(new['Category2'].value_counts())
		cat3 = pd.DataFrame(new['Category3'].value_counts())
		cat1.reset_index(inplace=True)
		cat2.reset_index(inplace=True)
		cat3.reset_index(inplace=True)
		cat1.columns = ['Category','Count'] #renaming columns
		cat2.columns = ['Category','Count']
		cat3.columns = ['Category','Count']
		cat1['Percent'] = (cat1['Count']/sum(cat1['Count']))*100 #Add column for percent
		cat2['Percent'] = (cat2['Count']/sum(cat2['Count']))*100
		cat3['Percent'] = (cat3['Count']/sum(cat3['Count']))*100
		# Write to	file for reading from R
		feather.write_dataframe(cat1, "cat1.feather")
		feather.write_dataframe(cat2, "cat2.feather")
		feather.write_dataframe(cat3, "cat3.feather")
		print("feather objects were saved")
	if args.countoption == 'taxa':
		#pull out columns with nitrogen fixing KOs
		Nif_KOs = new2[new2['accession_x'].str.match("K02588|K02586|K02591")] #nifHDK
		Nif_KOsT = pd.DataFrame(Nif_KOs['Taxa'].value_counts())
		Nif_KOsT.columns = ['Number of Gene Calls with NifHDK KOs']
		Nif_KOsT.index.name = 'Taxa'
		Nif_KOsT.reset_index(inplace=True)
		print("The Taxa are:")
		print(Nif_KOsT)
		Nif_KOs = pd.DataFrame(Nif_KOs['Lower_taxa'].value_counts())
		Nif_KOs.columns = ['Number of Gene Calls with NifHDK KOs']
		Nif_KOs.index.name = 'Lower_Taxa'
		Nif_KOs.reset_index(inplace=True)
		print("The Lower Taxa are:")
		print(Nif_KOs)
	else:
		print("You need to pick either taxa or cat for counts")
