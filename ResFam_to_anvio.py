#importing packages
import glob,os
import pandas as pd
import xlrd
#making txt file for gene.txt file to run anvi-run-hmm
#setting working directory, change folder to the location of your tblout files
#os.chdir("C:/Users/jvhagey/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/Foam_HMMER/")
os.chdir("C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/ResFams")
#Start at row 3 for the header
df = pd.read_excel("180102_resfams_metadata_updated_v1.2.2.xlsx", usecols='A:B,F', header=3)
#Move columns around and rename them
df = df[['Resfam Family Name', 'ResfamID', 'HMM Source']]
df.columns = ["gene","accession","hmmsource"]
#write DataFrame to tab separated file (.csv) with file name and TIGRFAM counts
df.to_csv('Resfam_anvi_genes.txt', sep='\t',index=False)

#even after this there are a few changes that need to be made to the gene names
