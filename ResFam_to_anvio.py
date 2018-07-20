#In shell terminal grep to get name and accession numbers of hmms
#grep -A1 "NAME" FOAM-hmm_rel1a.hmm > Acc_num.txt
#Use zgrep for .gz file
#importing packages
import pandas as pd
import re
import glob,os
import csv
#making txt file for gene.txt file to run anvi-run-hmm
os.chdir("C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/ResFams")
#file = "C:/Users/Jill/Desktop/Acc_num.txt"
#lines = open(file).readlines() #read linies of file
with open("Acc_num.txt", "r+") as f: #open file for reading and writing
    filedata = f.read()
    #Removing lines starting with --
    filedata = re.sub(r'--\n', '', filedata, flags=re.MULTILINE)
    #replacing two or three spaces with only one
    filedata = re.sub(r'   ', ' ', filedata, flags=re.MULTILINE)
    filedata = re.sub(r'  ', ' ', filedata, flags=re.MULTILINE)
    #turns out that resfam hmms have spaces between NAME and xxx and others have spaces
    filedata = re.sub(r' ', '\t', filedata, flags=re.MULTILINE)
    #One hmm has \t\s between NAME and name so fixing that next
    filedata = re.sub(r'\t\t', '\t ', filedata, flags=re.MULTILINE)
    filedata = re.sub(r' ', '', filedata, flags=re.MULTILINE)
with open("Acc_num.txt", 'w') as f2: #writing new header for edited file removing blank lines at beginning and end
    f2.writelines(filedata)
df = pd.read_csv("Acc_num.txt", sep='\t', header=None)
df.columns = ['Columns','Rows'] #rename columns
df = df.pivot(columns = "Columns", values = "Rows")
#remove none values and move up the cells
df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')
df = df[['NAME', 'ACC']]
df["hmmsource"] = "ResFams"
df.columns = ["gene","accession","hmmsource"] #rename columns
#write DataFrame to tab separated file (.csv)
df.to_csv('ResFam_gene_2.txt', sep='\t',index=False)
