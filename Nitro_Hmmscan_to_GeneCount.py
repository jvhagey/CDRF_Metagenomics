#Script takes tblout output file from hmmscan and converts them to a table that displays gene counts that can be transfered into R for analysis/graphing
#importing packages
import glob,os
import pandas as pd
import re
import csv
#setting working directory
os.chdir("C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/PyTest")
#saving new header to be used later in rewriting file
new_header = ("target_name accession query_name E-value score bias E-value score bias exp reg clu ov env dom rep inc description_of_target")
#files with extention that will be looped through code
hmmer_files = glob.glob("*.Nitro.tblout")
csvfile = "C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/PyTest/GeneCountOutput_Nitro.csv"
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
    filedata = filedata.replace(' [', '_')
    filedata = filedata.replace('] ', '_')
    filedata = filedata.replace('-[', '_')
    filedata = filedata.replace('"', '')
    #Looks for pattern of any letter except "T" with a single space then another letter except "T" and replaces space with _
    filedata = re.sub(r'([A-SU-Za-z])\s([A-SU-Za-z])', r'\1_\2',filedata)
    with open(file, 'w') as f2: #writing new header for edited file removing blank lines at beginning and end
        f2.seek(0)
        f2.write(new_header)
        f2.write('\n')
        f2.writelines(filedata[1:-1])
    #read in only columns I want to keep for analysis
    df = pd.read_csv(file, delim_whitespace=True, header=0, usecols=["target_name","E-value","description_of_target"])
    df_less = pd.DataFrame(df[df['E-value'] > 1.0E-15]) #Make dataframe with E-values less than 1x10^-15
    counts_DF = df_less.target_name.value_counts().reset_index().rename(columns={'index': 'Target_Name', 0: 'Gene_Count'})
    counts_DF['Sample_Name'] = open(file,'r').name #getting name of file currently looping and add to column
    counts_DF.columns = ['target_name','Gene_Count','Sample_Name']
    RecA_count = counts_DF.loc[counts_DF["target_name"].str.contains('TIGR02012'), "Gene_Count"].values[0]
    counts_DF['RecA_Count'] = RecA_count
    counts_DF_names = pd.merge(counts_DF, df_less[['target_name','description_of_target']], on='target_name').drop_duplicates()
    appended_counts.append(counts_DF_names) #store DataFrame in list
appended_counts = pd.concat(appended_counts) #combing output of for loop into one dataframe
appended_counts['Sample_Name'] = appended_counts['Sample_Name'].str.replace('_S(.*?).Nitro.tblout', '') #editing "Sample_Names" to merge libray size
appended_counts['Sample_Name'] = appended_counts['Sample_Name'].str.replace('-', '_')
#Run in server on output of flash2 to get library size used in hmmscan
#grep "Combined pairs" *.out > Merged_Stats.txt
Merged_Stats = pd.read_csv("Merge_Stats.txt", header=0, usecols=["Sample","Combined pairs"])
Merged_Stats_df = pd.DataFrame(data=Merged_Stats)
Merged_Stats_df.columns = ['Sample_Name', 'Combined_Pairs']
Norm_Counts_df = pd.merge(appended_counts, Merged_Stats_df, on='Sample_Name')
#Normalizing gene count to reads per million (RPM)
Norm_Counts_df["RPM"] = (Norm_Counts_df["Gene_Count"]/Norm_Counts_df["Combined_Pairs"])*1000000
#Normaling by RecA count
Norm_Counts_df["RecA_Norm"] = Norm_Counts_df["Gene_Count"]/Norm_Counts_df["RecA_Count"]
#Shortening description to gene name
Norm_Counts_df['Gene_Name'] = pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifH", case=False), "nifH",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifL", case=False), "nifL",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifA", case=False), "nifA",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifD", case=False), "nifD",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifX", case=False), "nifX",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifM", case=False), "nifM",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifB", case=False), "nifB",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifE", case=False), "nifE",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifV", case=False), "nifV",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifK", case=False), "nifK",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifU", case=False), "nifU",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("RecA", case=False), "recA",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("NifN", case=False), "nifN",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("VnfD", case=False), "vnfD",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("VnfG", case=False), "vnfG",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("VnfK", case=False), "vnfK",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("AnfK", case=False), "anfK",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("AnfD", case=False), "anfD",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("AnfG", case=False), "anfG",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("fdxN", case=False), "fdxN",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("alt_nitrog_alph", case=False), "nitrogenase alpha chain",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("dinitro_DRAG", case=False), "dinitrogen reductase",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("TIGR02935", case=False), "probable nitrogen fixation protein",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("dinitro_DRAG", case=False), "dinitrogen reductase",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("N2-ase-Ialpha", case=False), "nitrogenase component I, alpha chain", "No_name")))))))))))))))))))))))))
Norm_Counts_df['Farm'] = pd.np.where(Norm_Counts_df.Sample_Name.str.contains("8"), "Farm_8",
                   pd.np.where(Norm_Counts_df.Sample_Name.str.contains("6"), "Farm_6",
                   pd.np.where(Norm_Counts_df.Sample_Name.str.contains("5"), "Farm_5",
                   pd.np.where(Norm_Counts_df.Sample_Name.str.contains("1"), "Farm_1", "No_name"))))
#write DataFrame to comma separated file (.csv) with file name and TIGRFAM counts
Norm_Counts_df.to_csv('GeneCountOutput_Nitro.csv', sep=',')
