#importing packages
import pandas as pd
#making txt file for gene.txt file to run anvi-run-hmm
df = pd.read_excel("C:/Users/jvhagey/Downloads/180102_resfams_metadata_updated_v1.2.2.xlsx", usecols='A:B,F', header=0)
df.columns = ["gene","accession","hmmsource"]
#write DataFrame to tab separated file (.csv) with file name and TIGRFAM counts
df.to_csv('C:/Users/jvhagey/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/Anvio/Resfam_anvi_genes.txt', sep='\t',index=False)
