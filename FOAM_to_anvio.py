#In shell terminal grep to get name and accession numbers of hmms
#grep -A1 "NAME" FOAM-hmm_rel1a.hmm > Acc_num.txt
#importing packages
import pandas as pd
import re
#making txt file for gene.txt file to run anvi-run-hmm
file = "C:/Users/Jill/Desktop/Acc_num.txt"
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
#write DataFrame to tab separated file (.csv)
df.to_csv('C:/Users/Jill/Desktop/gene_FOAM.txt', sep='\t',index=False)
