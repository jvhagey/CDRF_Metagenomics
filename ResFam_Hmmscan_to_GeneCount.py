#importing packages
import glob,os
import pandas as pd
import re
import csv
#setting working directory
os.chdir("C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/PyTest")
#saving new header to be used later in rewriting file
new_header = ("target_name ResfamID query_name E-value score bias E-value score bias exp reg clu ov env dom rep inc description_of_target ARO")
#files with extention that will be looped through code
hmmer_files = glob.glob("*.ResFam.tblout")
csvfile = "C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/PyTest/GeneCountOutput_Resfam.csv"
ResFam_Meta = "C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/ResFams/resfams_metadata_v1.2.2.csv"
Nitro_Counts = "C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/PyTest/GeneCountOutput_Nitro.csv"
appended_counts = [] #making blank dataframe to store counts of each gene from output of for loop
#this loops through *.tblout files from HMMER and removes comment lines so that data analysis can continue down stream
for file in hmmer_files:
    lines = open(file).readlines() #read linies of file
    with open(file, "r+") as f: #open file for reading and writing
        f.seek(0)   #move pointer to top of file
        f.write('\n') #write a new line
        f.writelines(lines[3:-11]) #write lines removing the first 3 and last 11
    with open(file, 'r') as f2: #Read in the new file that was just written
        filedata = f2.read()
    #Replacing special characters in description_of_target and leaving only one space between coloumns. This are specific to Nitrofiles currently
    filedata = re.sub("\s\s+" , " ", filedata)
    filedata = filedata.replace(' -', '')
    filedata = filedata.replace(', ', '_')
    filedata = filedata.replace(': ', '_')
    filedata = filedata.replace(') ', '_')
    filedata = filedata.replace(' (', '_')
    filedata = filedata.replace(')', '')
    filedata = filedata.replace('"', '')
    filedata = filedata.replace('/', '_')
    filedata = filedata.replace('[', '')
    filedata = filedata.replace(']', '')
    #removes space between TIGRFAM number and PF numbers and beginning of description_of_target
    filedata = re.sub(r'(TIGR\d+)\s([A-Za-z])', r'\1_\2',filedata)
    filedata = re.sub(r'(PF\d+.\d+)\s([A-Za-z])', r'\1_\2',filedata)
    filedata = re.sub(r'([A-EG-QS-Za-z]) ([0-9])', r'\1_\2',filedata)
    #removing spaces between words that above code did not encompass
    filedata = filedata.replace('ABC1 family','ABC1_family')
    filedata = filedata.replace('PC1 beta','PC1_beta')
    filedata = filedata.replace('L1 beta-lactamase', 'L1_beta-lactamase')
    filedata = filedata.replace('class I and','class I_and')
    filedata = filedata.replace('ABC-2 type','ABC-2_type')
    filedata = filedata.replace('APH3''','APH3')
    filedata = filedata.replace('drug:H+ antiporter-2_14 Spanner','drug_H+_antiporter-2_14_Spanner')
    #Looks for pattern of any letter with a single space then another letter and replaces space with _
    filedata = re.sub(r'([A-Za-z]) ([A-Za-z])', r'\1_\2',filedata)
    #Putting spaces back where I want them for the second and last column
    filedata = re.sub(r'([A-Za-z])_(ARO:\d+)', r'\1 \2',filedata)
    filedata = re.sub(r'([A-Za-z0-9])_(RF\d+)', r'\1 \2',filedata)
    #Cleaning up a few more lines
    filedata = filedata.replace('TIGR01205_D_ala_D alaTIGR_D-alanine--D-alanine_ligase','TIGR01205_D_ala_D_alaTIGR_D-alanine--D-alanine_ligase')
    filedata = filedata.replace('Class_A beta-lactamase','Class_A_beta-lactamase')
    filedata = filedata.replace('Class_B beta-lactamase','Class_B_beta-lactamase')
    filedata = filedata.replace('Class_C beta-lactamases','Class_C_beta-lactamases')
    filedata = filedata.replace('Class_D beta-lactamases','Class_D_beta-lactamases')
    with open(file, 'w') as f2: #writing new header for edited file removing blank lines at beginning and end
        f2.seek(0)
        f2.write(new_header)
        f2.write('\n')
        f2.writelines(filedata[1:-1])
    #read in only columns I want to keep for analysis
    df = pd.read_csv(file, delim_whitespace=True, header=0, usecols=["target_name","ResfamID","E-value","description_of_target","ARO"])
    df_less = pd.DataFrame(df[df['E-value'] > 1.0E-15]) #Make dataframe with E-values less than 1x10^-15
    counts_DF = df_less.target_name.value_counts().reset_index().rename(columns={'index': 'target_name', 0: 'Gene_Count'})#Getting counts of Resfam_Family_Name
    counts_Res_DF = df_less.ResfamID.value_counts().reset_index().rename(columns={'index': 'ResfamID', 0: 'Gene_Count'})#Getting counts of Resfam_Family_Name
    #print(counts_DF.shape[0]) #gives number of row count
    #print(counts_Res_DF.shape[0]) #gives number of row count
    counts_DF['Sample_Name'] = open(file,'r').name #getting name of file currently looping and add to column
    counts_DF.columns = ['target_name','Gene_Count','Sample_Name'] #rename columns
    counts_DF_names = pd.merge(counts_DF, df_less[['target_name','ResfamID','description_of_target']], on='target_name').drop_duplicates()
    appended_counts.append(counts_DF_names) #store DataFrame in list
appended_counts = pd.concat(appended_counts) #combing output of for loop into one dataframe
appended_counts['Sample_Name'] = appended_counts['Sample_Name'].str.replace('_S(.*?).ResFam.tblout', '') #editing "Sample_Names" to merge libray size
appended_counts['Sample_Name'] = appended_counts['Sample_Name'].str.replace('-', '_')
#Importing gene count from Nitro hmms to get RecA counts for normalization
Nitro_count = pd.read_csv(Nitro_Counts, sep=',', header=0, usecols=["Gene_Count","Sample_Name","Gene_Name"])
RecA_count = pd.DataFrame(Nitro_count.loc[Nitro_count["Gene_Name"].str.contains('recA')])
appended_counts2 = pd.merge(appended_counts, RecA_count[["Gene_Count","Sample_Name"]], on='Sample_Name')#Adding new column with recA counts for each sample
appended_counts2.columns = ['Resfam_Family_Name','Gene_Count','Sample_Name','ResfamID','description_of_target','RecA_Count']
#Run in server on output of flash2 to get library size used in hmmscan
#grep "Combined pairs" *.out > Merged_Stats.txt
Merged_Stats = pd.read_csv("Merge_Stats.txt", header=0, usecols=["Sample","Combined pairs"])
Merged_Stats_df = pd.DataFrame(data=Merged_Stats)
Merged_Stats_df.columns = ['Sample_Name', 'Combined_Pairs']
Norm_Counts_df = pd.merge(appended_counts2, Merged_Stats_df, on='Sample_Name')
#Normalizing gene count to reads per million (RPM)
Norm_Counts_df["RPM"] = (Norm_Counts_df["Gene_Count"]/Norm_Counts_df["Combined_Pairs"])*1000000
#Normaling by RecA count
Norm_Counts_df["RecA_Norm"] = Norm_Counts_df["Gene_Count"]/Norm_Counts_df["RecA_Count"]
#Adding Metadata (type of ) to dataframe on ResFams
ResFam_Meta_df = pd.read_csv(ResFam_Meta, sep=',', header=0, usecols=["Resfam Family Name","Antibiotic Classification (Resfam Only)","Mechanism Classification","Beta-Lactamase Ambler Class","Tetracycline Resistance (Resfam Only)"])
ResFam_Meta_df.columns = ["Resfam_Family_Name","Antibiotic_Classification","Mechanism_Classification","Beta-Lactamase_Ambler_Class","Tetracycline_Resistance"]
Counts_Meta_df = pd.merge(Norm_Counts_df, ResFam_Meta_df, on='Resfam_Family_Name')
Counts_Meta_df['Farm'] = pd.np.where(Counts_Meta_df.Sample_Name.str.contains("8"), "Farm_8",
                   pd.np.where(Counts_Meta_df.Sample_Name.str.contains("6"), "Farm_6",
                   pd.np.where(Counts_Meta_df.Sample_Name.str.contains("5"), "Farm_5",
                   pd.np.where(Counts_Meta_df.Sample_Name.str.contains("1"), "Farm_1", "No_name"))))
#write DataFrame to comma separated file (.csv) with file name and TIGRFAM counts
Counts_Meta_df.to_csv('GeneCountOutput_Resfam.csv', sep=',')
