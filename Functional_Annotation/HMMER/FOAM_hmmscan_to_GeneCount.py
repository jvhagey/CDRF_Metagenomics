#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019
## Given a hmmr .tblout file from searching against FOAM HMMs, script will give a .feather file to read into R
#Script takes tblout output file from hmmscan and converts them to a table that displays gene counts that can be transfered into R for analysis/graphing

#importing packages
import glob,os
import pandas as pd
import re
import csv
import feather
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="FOAM_hmmscan_to_GeneCount.py",description="Parses Hmmscan .tblout files into an R readable .feather file")
    parser.add_argument("-s", dest="stats", action="store", required=True, help='Give name of Stats_File.txt file. Should be in the same folder as hmmscan files. This is a tab separated file with the sample name and read/contig count (required)')
    #parser.add_argument("-nk", dest="nitroKo", action="store", default=None, required=True, help="Give name of Nitrogen_cycle_KO.tsv file. (required)")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were .tblout files are found (required)")
    #parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None, required=True, help="Output directory (required)")
    return parser.parse_args()

if __name__ == '__main__':
    #Parse command-line
    args = parse_cmdline()
    #Make output directory
    #os.makedirs(args.outdirname)
    #saving new header to be used later in rewriting file
    new_header = ("target_name accession query_name E-value score bias E-value score bias exp reg clu ov env dom rep inc description_of_target")
    #Calling in files
    os.chdir(args.indirname) #change directory
    hmmer_files = glob.glob("*.FOAM.tblout") #files with extention that will be looped through code
#csvfile = "C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/PyTest/GeneCountOutput_Nitro.csv"
    appended_counts = [] #making blank dataframe to store counts of each gene from output of for loop
    #this loops through *.tblout files from HMMER and removes comment lines so that data analysis can continue down stream
    for file in hmmer_files:
        with open(file, 'r') as f2: #Read in the new file that was just written
            filedata = f2.read()
            filedata = re.sub(r'(?m)^\#.*\n?' , "", filedata)
        #Replacing special characters in description_of_target and leaving only one space between coloumns.
        filedata = re.sub("\s\s+" , " ", filedata)
        filedata = filedata.replace(' -', '')
        filedata = filedata.replace(', ', '_')
        #removes space between KO number and pathway numbers and beginning of description_of_target
        filedata = re.sub(r'(K\d+) (\d+.)', r'\1_\2',filedata)
        with open(file, 'w') as f2: #writing new header for edited file removing blank lines at beginning and end
            f2.seek(0)
            f2.write(new_header)
            f2.write('\n')
            f2.writelines(filedata)
        #read in only columns I want to keep for analysis
        df = pd.read_csv(file, delim_whitespace=True, header=0, usecols=["target_name","E-value","description_of_target"])
        df_less = pd.DataFrame(df[df['E-value'] > 1.0E-15]) #Make dataframe with E-values less than 1x10^-15
        counts_DF = df_less.target_name.value_counts().reset_index().rename(columns={'index': 'Target_Name', 0: 'Gene_Count'})
        counts_DF['Sample_Name'] = open(file,'r').name #getting name of file currently looping and add to column
        counts_DF.columns = ['target_name','Gene_Count','Sample_Name']
        counts_DF_names = pd.merge(counts_DF, df_less[['target_name','description_of_target']], on='target_name').drop_duplicates()
        appended_counts.append(counts_DF_names) #store DataFrame in list
appended_counts = pd.concat(appended_counts) #combing output of for loop into one dataframe
appended_counts['Sample_Name'] = appended_counts['Sample_Name'].str.replace('.FOAM.tblout', '') #editing "Sample_Names" to merge libray size
#Run in server on output of megahit_tools to get library size used in hmmscan
#grep "[STAT]" *.out > Merged_Stats.txt
Merged_Stats = pd.read_csv(args.stats, header=0, usecols=["Sample_Name","Contigs"])
Merged_Stats_df = pd.DataFrame(data=Merged_Stats)
Norm_Counts_df = pd.merge(appended_counts, Merged_Stats_df, on='Sample_Name')
#Normalizing gene count to reads per million (RPM)
Norm_Counts_df["RPM"] = (Norm_Counts_df["Gene_Count"]/Norm_Counts_df["Contigs"])*1000000
#Numbers at beginning of file are Month-Farm-Sample
Norm_Counts_df['Farm'] = pd.np.where(Norm_Counts_df.Sample_Name.str.contains("1"), "Farm_1", #Adding column that will have the farm name
                   pd.np.where(Norm_Counts_df.Sample_Name.str.contains("6"), "Farm_6",
                   pd.np.where(Norm_Counts_df.Sample_Name.str.contains("2"), "Farm_2",
                   pd.np.where(Norm_Counts_df.Sample_Name.str.contains("1"), "Farm_1", "No_name"))))
Norm_Counts_df['Month'] = pd.np.where(Norm_Counts_df.Sample_Name.str.contains("12"), "December", #Adding column that will have the Month
                pd.np.where(Norm_Counts_df.Sample_Name.str.contains("10"), "October",
                pd.np.where(Norm_Counts_df.Sample_Name.str.contains("4"), "April",
                pd.np.where(Norm_Counts_df.Sample_Name.str.contains("7"), "July", "No_name"))))
                #Download FOAM-onto_rel1.tsv from FOAM database and run the following line on linux to get nitrogen cycle related KOs
                #grep "11_Nitrogen cycle" FOAM-onto_rel1.tsv > Nitrogen_cycle_KO.tsv
Norm_Counts_df.to_csv('Norm_Counts_df.csv', sep=',',index=False)
#Shortening description to gene name
Norm_Counts_df["Function"] = pd.np.where(Norm_Counts_df.description_of_target.str.contains("K00370|K00371|K00372|K00373|K00374|K00345|K00346|K00347|K08361|K02567|K02568|K15878|K15879", case=False), "Denitrification: Nitrate to nitrite",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K10534|K00360|K15875", case=False), "Assimilatory nitrate reduction: Nitrate to nitrite",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K00368|K15864", case=False), "Denitrification: Nitrite to nitrite oxide",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K04561|K02305", case=False), "Denitrification: Nitrite oxide to nitrous oxide",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K00376", case=False), "Denitrification: Nitrous oxide to nitrogen",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K03385|K04013|K04014|K04015", case=False), "Dissimilatory nitrate reduction: Nitrate to ammonia",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K00362|K00366|K00363", case=False), "Assimilatory nitrate reduction: Nitrate to ammonia",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K05601|K15864", case=False), "Hydroxylamine reduction",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K10944|K10945|K10946", case=False), "Ammonia to hydroxylamine",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K10535", case=False), "Hydroxylamine to nitrite",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K02588|K02591|K02586|K00531|K02595", case=False), "Nitrogen Fixation",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K00260|K00261|K00262|K01915|K15371|K00265|K00266|K00264|K05597|K01425|K05597|K01953|K00284|K00264", case=False), "Ammonia assimilation",
                   pd.np.where(Norm_Counts_df.description_of_target.str.contains("K15905|K15906", case=False), "Nitrification: Nitrite to nitrate","No_name")))))))))))))
#write DataFrame to comma separated file (.csv) with file name and FOAM hmm counts
Norm_Counts_df.to_csv('GeneCountOutput_Nitro.csv', sep=',',index=False)
#Write to file for reading from R
feather.write_dataframe(Norm_Counts_df, "Hmm_hit_counts.feather")
