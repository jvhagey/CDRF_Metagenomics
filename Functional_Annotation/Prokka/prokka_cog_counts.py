#importing packages
import feather
import glob,os
import pandas as pd
import re
import csv
#setting working directory
os.chdir("/mnt/c/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/COGs")
#making txt file
file_cat = "cog2cat23814_16-nov-2018.csv"
file_cogs = "Standard.cog"
#making dataframes
df_cat = pd.read_csv(file_cat, sep=',', header=0)
df_cogs = pd.read_csv(file_cogs, sep='\t', header=None)
df_cogs.columns = ['Gene_Call','COG ID'] 
df = pd.merge(df_cat, df_cogs, on='COG ID', how="right")
#df.to_csv('cogs_and_cats.csv', sep=',' ,index=False) # writing to a .csv file
#getting frequence counts
counts = pd.DataFrame(df['Catergory'].value_counts())
#Write to  file for reading from R
feather.write_dataframe(counts, "counts.feather")