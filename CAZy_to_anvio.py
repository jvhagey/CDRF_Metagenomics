#importing packages
import pandas as pd
#making txt file for gene.txt file to run anvi-run-hmm
file = "C:/Users/jvhagey/Desktop/CAZyDB-ec-info.txt"
fields = ["#family","genbank"]
#make pandas DataFrame with only the first two columns of data
df = pd.read_csv(file, sep='\t', lineterminator='\n', usecols=fields, header=0)
df = df[["#family","genbank"]]
#Add in 3rd column
df['hmmsource'] = "CAZyDB"
df.columns = ["gene","accession","hmmsource"]
#write DataFrame to tab separated file (.csv)
df.to_csv('C:/Users/jvhagey/Desktop/gene.txt', sep='\t',index=False)
