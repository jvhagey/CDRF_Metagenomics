#importing packages
import feather
import pandas as pd
import re
import glob,os
import csv
#setting working directory
os.chdir("C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/GhostKOALA")
#making txt file
file_ortho = "KeggOrthology_Table1.txt"
file_userdf = "user_ko_definition_all.txt"
file_user = "user_ko_all.txt"
#making dataframes
df_ortho = pd.read_csv(file_ortho, sep=',', header=0)
df_userdf = pd.read_csv(file_userdf, sep='\t', header=None)
df_userdf.columns = ['Query','accession','Definition','Score','Second best','Second Score'] #renaming columns
df_userdf = df_userdf.drop(columns=['Definition','Score','Second best','Second Score']) #droping columns
new = pd.merge(df_userdf, df_ortho, on='accession')
cat1 = pd.DataFrame(new['Category1'].value_counts())
cat2 = pd.DataFrame(new['Category2'].value_counts())
cat3 = pd.DataFrame(new['Category3'].value_counts())
cat1.reset_index(inplace=True)
cat2.reset_index(inplace=True)
cat3.reset_index(inplace=True)
cat1.columns = ['Category','Count'] #renaming columns
cat2.columns = ['Category','Count']
cat3.columns = ['Category','Count']
cat1['Percent'] =  (cat1['Count']/sum(cat1['Count']))*100 #Add column for percent 
cat2['Percent'] =  (cat2['Count']/sum(cat2['Count']))*100
cat3['Percent'] =  (cat3['Count']/sum(cat3['Count']))*100
# Write to  file for reading from R
feather.write_dataframe(cat1, "cat1.feather")
feather.write_dataframe(cat2, "cat2.feather")
feather.write_dataframe(cat3, "cat3.feather")
