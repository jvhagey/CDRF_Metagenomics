#importing packages
import feather
import glob,os
import pandas as pd
from pandas import DataFrame
import re
import csv
#setting working directory
os.chdir("/mnt/c/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/COGs")
#making dataframes from csv files
df_cat = pd.read_csv('Cog_Cats.txt', sep='\t', header=0)
df_let = pd.read_csv("cognames2003-2014.tab", sep='\t', header=0, encoding='windows-1252')
df_cogs = pd.read_csv("Standard.cog", sep='\t', header=None)
#changing column names
df_cogs.columns = ['Gene_Call','COG ID']
df_let.columns = ['COG_ID','Letter','Name']
#adding addtional rows for cases where there are two letters in letter column
##first we separate multiple capitial letters in the "Letters" columns by commas by applying a function across the dataframe
def split_it(x):
    return re.sub(r'([A-Z])([A-Z])', r'\1,\2',x)
    return re.sub(r'([A-Z])([A-Z])([A-Z])', r'\1,\2,\3',x)
df_let['Letter'] = df_let['Letter'].apply(split_it)
new_let = DataFrame(df_let.Letter.str.split(',').tolist(), index=([df_let.Name,df_let.COG_ID])).stack()
new_let = new_let.reset_index()[[0, 'COG_ID', 'Name']] # "Letter" variable is currently labeled 0
new_let.columns = ['Letter','COG ID','Name'] # renaming Letter
#Now we can move one to combining dataframe
df = pd.merge(df_cat, new_let, on='Letter', how="right")
df.to_csv('cogs_and_cats.csv', sep=',' ,index=False) # writing to a .csv file
df_cogs_new = pd.merge(df, df_cogs, on='COG ID', how="right")
df_cogs_new.to_csv('genecalls_to_cats.csv', sep=',' ,index=False) # writing to a .csv file
#Note that df_cogs_new will have a different shape than df_cogs since some cogs have multiple categories assigned to them
#getting frequency counts
counts = pd.DataFrame(df['Catergory'].value_counts())
counts.reset_index(inplace=True)
counts.columns = ['Catergory','Counts']
counts['Percent'] =  (counts['Counts']/sum(counts['Counts']))*100
#Write to  file for reading from R
feather.write_dataframe(counts, "Prokka_cog_counts.feather")
