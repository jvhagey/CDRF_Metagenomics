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
    parser.add_argument("-s", dest="stats", action="store", required=True, help='Give name of Stats_File.txt ("Merge_Stats.txt") file. Should be in the same folder as hmmscan files. This is a tab separated file with the sample name and read/contig count (required)')
    parser.add_argument("--recA", dest="recA", action="store", required=True, help='Give name of file with RecA count data ("GeneCountOutput_Nitro.csv"). Should be in the same folder as hmmscan files. This is a tab separated file with the sample name and read/contig count (required)')
    parser.add_argument("-r", dest="results", action="store", required=True, help='Give name of Results_File.txt with genome size data ("Merged_Results.txt") file. Should be in the same folder as hmmscan files. This is a tab separated file with the sample name and read/contig count (required)')
    #parser.add_argument("-nk", dest="nitroKo", action="store", default=None, required=True, help="Give name of Nitrogen_cycle_KO.tsv file. (required)")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were .tblout files are found (required)")
    #parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None, required=True, help="Output directory (required)")
    return parser.parse_args()

def clean_up(outfile):
    os.remove(outfile)

def main():
    #Parse command-line
    args = parse_cmdline()
    #Foam_Meta = "C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/dairy sequencing/Metagenomics/Foam_HMMER/FOAM-onto_rel1.tsv"
    #Nitro_Counts = "C:/Users/Jill/OneDrive - UC Davis/Documents/collaboration/Metagenomics_2016/PyTest/GeneCountOutput_Nitro.csv"
    #Merge_Stats_File= "C:/Users/jvhagey/OneDrive - UC Davis/Documents/collaboration/Metagenomics_2016/HMMER_Copy/Merge_Stats.txt"
    #saving new header to be used later in rewriting file
    new_header = ("target_name accession query_name E-value score bias E-value score bias exp reg clu ov env dom rep inc description_of_target")
    #Calling in files
    os.chdir(args.indirname) #change directory
    hmmer_files = glob.glob("*.Foam.tblout") #files with extention that will be looped through code
    appended_counts = [] #making blank dataframe to store counts of each gene from output of for loop
    #this loops through *.tblout files from HMMER and removes comment lines so that data analysis can continue down stream
    for file in hmmer_files:
        with open(file, 'r') as infile, open(file + "_temp.txt", 'w') as outfile: #Read in the new file that was just written
            temp = infile.read().replace(' -', '')
            temp = re.sub("\s\s+" , " ", temp)
            temp = temp.replace(', ', '_')#Replacing special characters in description_of_target and leaving only one space between coloumns.
            #temp = re.sub(r'(?m)^\#.*\n?' , "", temp)
            temp = re.sub(r'(K\d+) (\d+.)', r'\1_\2',temp)#removes space between KO number and pathway numbers and beginning of description_of_target
            outfile.write(temp)
        with open(file +"_temp.txt", 'r') as infile, open(file+"_clean.txt", "w") as outfile : #writing new header for edited file removing blank lines at beginning and end
            outfile.write(new_header)
            outfile.write('\n')
            for line in infile:
                if not line.startswith('#'):
                    outfile.write(line)
        clean_up(file+"_temp.txt")
        #read in only columns I want to keep for analysis
        df = pd.read_csv(file + "_clean.txt", delim_whitespace=True, header=0, usecols=["target_name","E-value","description_of_target"])
        df_less = pd.DataFrame(df[df['E-value'] > 1.0E-15]) #Make dataframe with E-values less than 1x10^-15
        counts_DF = df_less.target_name.value_counts().reset_index().rename(columns={'index': 'Target_Name', 0: 'Gene_Count'})
        counts_DF['Sample_Name'] = open(file,'r').name #getting name of file currently looping and add to column
        counts_DF.columns = ['target_name','Gene_Count','Sample_Name']
        counts_DF_names = pd.merge(counts_DF, df_less[['target_name']], on='target_name').drop_duplicates()
        appended_counts.append(counts_DF_names) #store DataFrame in list
        print("\n"+'\033[93m' + "Finished cleaning " + open(file,'r').name +". New file created: " + open(file,'r').name +"_cleaned.txt")
    appended_counts = pd.concat(appended_counts) #combing output of for loop into one dataframe
    appended_counts['Sample_Name'] = appended_counts['Sample_Name'].str.replace('_S(.*?).Foam.tblout', '') #editing "Sample_Names" to merge libray size
    appended_counts['Sample_Name'] = appended_counts['Sample_Name'].str.replace('-', '_') #editing "Sample_Names" to merge libray size
    Nitro_count = pd.read_csv(args.recA, sep=',', header=0, usecols=["Gene_Count","Sample_Name","Gene_Name"])
    RecA_count = pd.DataFrame(Nitro_count.loc[Nitro_count["Gene_Name"].str.contains('recA')])
    appended_counts2 = pd.merge(appended_counts, RecA_count[["Gene_Count","Sample_Name"]], on='Sample_Name') #Adding new column with recA counts for each sample
    appended_counts2.columns = ['KO_Name','Gene_Count','Sample_Name','RecA_Count']  #Run in server on output of megahit_tools to get library size used in hmmscan
    #grep "[STAT]" *.out > Merged_Stats.txt
    Merged_Stats = pd.read_csv(args.stats, header=0, usecols=["Sample","Combined pairs"])
    Merged_Stats_df = pd.DataFrame(data=Merged_Stats)
    Merged_Stats_df.columns = ["Sample_Name","Combined_Pairs"]
    Norm_Counts_df = pd.merge(appended_counts2, Merged_Stats_df, on='Sample_Name')
    #Normalizing gene count to reads per million (RPM)
    Norm_Counts_df["RPM"] = (Norm_Counts_df["Gene_Count"]/Norm_Counts_df["Combined_Pairs"])*1000000
    ##Normaling by RecA count
    Norm_Counts_df["RecA_Norm"] = Norm_Counts_df["Gene_Count"]/Norm_Counts_df["RecA_Count"]
    #Normalizing by genome size
    Genome_size = pd.read_csv(args.results, header=0, sep='\t',usecols=["Sample_Name","average_genome_size"])
    Genome_size["Sample_Name"] = Genome_size["Sample_Name"].str.replace('_S\d+', '')
    Genome_size['Sample_Name'] = Genome_size['Sample_Name'].str.replace('-', '_') #editing "Sample_Names" to merge libray size
    Norm_Counts_df = pd.merge(Norm_Counts_df, Genome_size, on='Sample_Name')
    #Numbers at beginning of file are Month-Farm-Sample
    Norm_Counts_df['Farm'] = pd.np.where(Norm_Counts_df.Sample_Name.str.startswith("1"), "Farm 1", #Adding column that will have the farm name
                pd.np.where(Norm_Counts_df.Sample_Name.str.startswith("6"), "Farm 6",
                pd.np.where(Norm_Counts_df.Sample_Name.str.startswith("5"), "Farm 5",
                pd.np.where(Norm_Counts_df.Sample_Name.str.startswith("8"), "Farm 8", "No name"))))
                #Download FOAM-onto_rel1.tsv from FOAM database and run the following line on linux to get nitrogen cycle related KOs
                #grep "11_Nitrogen cycle" FOAM-onto_rel1.tsv > Nitrogen_cycle_KO.tsv
    Norm_Counts_df["Function"] = pd.np.where(Norm_Counts_df.KO_Name.str.contains("K00370|K00371|K00372|K00373|K00374|K00345|K00346|K00347|K08361|K02567|K02568|K15878|K15879", case=False), "Denitrification: Nitrate to nitrite",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K10534|K00360|K15875", case=False), "Assimilatory nitrate reduction: Nitrate to nitrite",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K00368|K15864", case=False), "Denitrification: Nitrite to nitrite oxide",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K04561|K02305", case=False), "Denitrification: Nitrite oxide to nitrous oxide",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K00376", case=False), "Denitrification: Nitrous oxide to nitrogen",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K03385|K04013|K04014|K04015", case=False), "Dissimilatory nitrate reduction: Nitrate to ammonia",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K00362|K00366|K00363", case=False), "Assimilatory nitrate reduction: Nitrate to ammonia",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K05601|K15864", case=False), "Hydroxylamine reduction",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K10944|K10945|K10946", case=False), "Ammonia to hydroxylamine",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K10535", case=False), "Hydroxylamine to nitrite",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K02588|K02591|K02586|K00531|K02595", case=False), "Nitrogen Fixation",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K00260|K00261|K00262|K01915|K15371|K00265|K00266|K00264|K05597|K01425|K05597|K01953|K00284|K00264", case=False), "Ammonia assimilation",
                   pd.np.where(Norm_Counts_df.KO_Name.str.contains("K15905|K15906", case=False), "Nitrification: Nitrite to nitrate","No_name")))))))))))))
    Norm_Counts_df['Function_L2'] = pd.np.where(Norm_Counts_df.KO_Name.str.contains("K00370|K00371|K00372|K10534|K00360|K15875|K08345|K00373|K00374|K08346|K08347|K08361|K02567|K02568|K15878|K15879|K15864|K04561|K02305|K00376|K03385|K04013|K04014|K04015|K00362|K00366|K00363|K05601|K15864", case=False), "Nitrogenated compounds Reduction",
                    pd.np.where(Norm_Counts_df.KO_Name.str.contains("K02588|K02591|K02586|K00531|K02595", case=False), "Nitrogen Fixation",
                    pd.np.where(Norm_Counts_df.KO_Name.str.contains("K10944|K10945|K10946|K10535|K15905|K15906", case=False), "Nitrification", "No_name")))
    Norm_Counts_df["Gene"] = pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00260_1.4.1.2", case=False), "gudB, rocG",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00261_1.4.1.3", case=False), "GLUD1_2, gdhA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00262_1.4.1.4", case=False), "gdhA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00264_1.4.1.13,1.4.1.14", case=False), "GLT1",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00265_1.4.1.13,1.4.1.14", case=False), "gltB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00265,KO:K00264_1.4.1.13,1.4.1.14", case=False), "GLT1, gltB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00265,KO:K00264,KO:K00284_1.4.1.13,1.4.1.14,1.4.7.1", case=False), "GLT1, gltB, GLU, gltS",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00265,KO:K00284_1.4.1.13,1.4.1.14,1.4.7.1", case=False), "GLT1, GLU, gltS",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00266_1.4.1.13,1.4.1.14", case=False), "gltD",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00266,KO:K00264_1.4.1.13,1.4.1.14", case=False), "gltD, GLT1",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00284_1.4.7.1", case=False), "GLU, gltS",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360_1.7.1.1", case=False), "nasB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360,KO:K00387_1.7.1.1,1.8.3.1", case=False), "nasB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360,KO:K00387,KO:K05301_1.7.1.1,1.8.3.1,1.8.2.1", case=False), "nasB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360,KO:K00507_1.7.1.1,1.14.19.1", case=False), "nasB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360,KO:K05301,KO:K00387_1.7.1.1,1.8.2.1,1.8.3.1", case=False), "nasB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362_1.7.1.4", case=False), "nirB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362,KO:K00363_1.7.1.4", case=False), "nirB, nirD",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362,KO:K00369_1.7.1.4,1.7.99.4", case=False), "nirB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362,KO:K00372_1.7.1.4,1.7.99.4", case=False), "nirB, nasA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362,KO:K00382_1.7.1.4,1.8.1.4", case=False), "nirB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00363_1.7.1.4", case=False), "nirD",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00363,KO:K00362_1.7.1.4", case=False), "nirD, nirB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00366_1.7.7.1", case=False), "nirA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00366,KO:K00392_1.7.7.1,1.8.7.1", case=False), "nirA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00366,KO:K00392,KO:K00381_1.7.7.1,1.8.7.1,1.8.1.2", case=False), "nirA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00366,KO:K01011_1.7.7.1,2.8.1.1,2.8.1.2", case=False), "nirA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368_1.7.2.1", case=False), "nirK",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K00405_1.7.2.1,1.9.3.1", case=False), "nirK",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K02275_1.7.2.1,1.9.3.1", case=False), "nirK",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K02277_1.7.2.1,1.9.3.1", case=False), "nirK",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K05601_1.7.2.1,1.7.99.1", case=False), "nirK, hcp",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K15864_1.7.2.1,1.7.99.1", case=False), "nirK, nirS",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K15864,KO:K05601_1.7.2.1,1.7.99.1", case=False), "nirK, nirS, hcp",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00370_1.7.99.4", case=False), "narG, narZ, nxrA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00370,KO:K08345_1.7.99.4", case=False), "narG, narZ, nxrA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00370,KO:K08345,KO:K15905_1.7.99.4", case=False), "narG, narZ, nxrA, nitrite oxidase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00371_1.7.99.4", case=False), "narH, narY, nxrB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00371,KO:K08346_1.7.99.4", case=False), "narH, narY, nxrB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00371,KO:K08346,KO:K15906_1.7.99.4", case=False), "narH, narY, nxrB, nitrite oxidase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00372_1.7.99.4", case=False), "nasA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00372,KO:K00369_1.7.99.4", case=False), "nasA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00372,KO:K02567_1.7.99.4", case=False), "nasA, napA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00372,KO:K02567,KO:K00369_1.7.99.4", case=False), "nasA, napA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00373", case=False), "narJ, narW",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00374_1.7.99.4", case=False), "narI, narV",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00374,KO:K08347_1.7.99.4", case=False), "narI, narV",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00376_1.7.2.4", case=False), "nosZ",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00376_1.7.99.6", case=False), "nosZ",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K01425_3.5.1.2", case=False), "glsA, GLS",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K01915_6.3.1.2", case=False), "glnA, GLUL",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K01953_6.3.5.4", case=False), "asnB, ASNS",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02305", case=False), "norC",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02305_1.7.99.7", case=False), "norC",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02567_1.7.99.4", case=False), "napA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02567,KO:K00372_1.7.99.4", case=False), "napA, nasA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02567,KO:K00372,KO:K00369_1.7.99.4", case=False), "napA, nasA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02568", case=False), "napB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02586_1.18.6.1", case=False), "nifD",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02588_1.18.6.1", case=False), "nifH",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02591_1.18.6.1", case=False), "nifK",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K03385_1.7.2.2", case=False), "nrfA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K04013", case=False), "nrfB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K04014", case=False), "nrfC",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K04015", case=False), "nrfD",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K04561_1.7.2.5", case=False), "norB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K05597_3.5.1.38", case=False), "aspQ, ansB, ansA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K05601_1.7.99.1", case=False), "hcp",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K05601,KO:K00368_1.7.99.1,1.7.2.1", case=False), "hcp, nirK",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K05601,KO:K15864,KO:K00368_1.7.99.1,1.7.2.1", case=False), "hcp, nirS, nirK",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K08347_1.7.99.4", case=False), "narI, narV",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K08361", case=False), "narJ, narW",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K10534_1.7.1.3", case=False), "NR",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K10535_1.7.3.4", case=False), "hao",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15864", case=False), "nirS",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15864_1.7.2.1,1.7.99.1", case=False), "nirS",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15864,KO:K00368_1.7.2.1,1.7.99.1", case=False), "nirS, nirK",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15864,KO:K00368,KO:K05601_1.7.2.1,1.7.99.1", case=False), "nirS, nirK, hcp",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15878", case=False), "narB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15878,KO:K03890", case=False), "narB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15906", case=False), "nitrite oxidase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15906,KO:K00371,KO:K08346_1.7.99.4", case=False), "nitrite oxidase, narH, narY, nxrB", "No_name"))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
    Norm_Counts_df["Description"] = pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00260_1.4.1.2", case=False), "glutamate dehydrogenase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00261_1.4.1.3", case=False), "glutamate dehydrogenase (NAD(P)+)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00262_1.4.1.4", case=False), "glutamate dehydrogenase (NADP+)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00264_1.4.1.13,1.4.1.14", case=False), "glutamate synthase (NADH)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00265_1.4.1.13,1.4.1.14", case=False), "glutamate synthase (NADPH) large chain",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00265,KO:K00264_1.4.1.13,1.4.1.14", case=False), "glutamate synthase (NADPH) large chain, glutamate synthase (NADH)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00265,KO:K00264,KO:K00284_1.4.1.13,1.4.1.14,1.4.7.1", case=False), "glutamate synthase (NADPH) large chain, glutamate synthase (NADH), glutamate synthase (ferredoxin)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00265,KO:K00284_1.4.1.13,1.4.1.14,1.4.7.1", case=False), "glutamate synthase (NADPH) large chain, glutamate synthase (ferredoxin)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00266_1.4.1.13,1.4.1.14", case=False), "glutamate synthase (NADPH) small chain",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00266,KO:K00264_1.4.1.13,1.4.1.14", case=False), "glutamate synthase (NADPH) small chain, glutamate synthase (NADH)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00284_1.4.7.1", case=False), "glutamate synthase (ferredoxin)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360_1.7.1.1", case=False), "assimilatory nitrate reductase electron transfer subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360,KO:K00387_1.7.1.1,1.8.3.1", case=False), "assimilatory nitrate reductase electron transfer subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360,KO:K00387,KO:K05301_1.7.1.1,1.8.3.1,1.8.2.1", case=False), "assimilatory nitrate reductase electron transfer subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360,KO:K00507_1.7.1.1,1.14.19.1", case=False), "assimilatory nitrate reductase electron transfer subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00360,KO:K05301,KO:K00387_1.7.1.1,1.8.2.1,1.8.3.1", case=False), "assimilatory nitrate reductase electron transfer subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362_1.7.1.4", case=False), "nitrite reductase (NADH) large subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362,KO:K00363_1.7.1.4", case=False), "nitrite reductase (NADH) large subunit, nitrite reductase (NADH) small subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362,KO:K00369_1.7.1.4,1.7.99.4", case=False), "nitrite reductase (NADH) large subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362,KO:K00372_1.7.1.4,1.7.99.4", case=False), "nitrite reductase (NADH) large subunit, assimilatory nitrate reductase catalytic subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00362,KO:K00382_1.7.1.4,1.8.1.4", case=False), "nitrite reductase (NADH) large subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00363_1.7.1.4", case=False), "nitrite reductase (NADH) small subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00363,KO:K00362_1.7.1.4", case=False), "nitrite reductase (NADH) small subunit, nitrite reductase (NADH) large subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00366_1.7.7.1", case=False), "ferredoxin-nitrite reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00366,KO:K00392_1.7.7.1,1.8.7.1", case=False), "ferredoxin-nitrite reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00366,KO:K00392,KO:K00381_1.7.7.1,1.8.7.1,1.8.1.2", case=False), "ferredoxin-nitrite reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00366,KO:K01011_1.7.7.1,2.8.1.1,2.8.1.2", case=False), "ferredoxin-nitrite reductase, ",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368_1.7.2.1", case=False), "nitrite reductase (NO-forming)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K00405_1.7.2.1,1.9.3.1", case=False), "nitrite reductase (NO-forming)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K02275_1.7.2.1,1.9.3.1", case=False), "nitrite reductase (NO-forming)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K02277_1.7.2.1,1.9.3.1", case=False), "nitrite reductase (NO-forming)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K05601_1.7.2.1,1.7.99.1", case=False), "nitrite reductase (NO-forming), hydroxylamine reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K15864_1.7.2.1,1.7.99.1", case=False), "nitrite reductase (NO-forming), nitrite reductase (NO-forming) / hydroxylamine reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00368,KO:K15864,KO:K05601_1.7.2.1,1.7.99.1", case=False), "nitrite reductase (NO-forming), nitrite reductase (NO-forming) / hydroxylamine reductase, hydroxylamine reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00370_1.7.99.4", case=False), "nitrate reductase / nitrite oxidoreductase, alpha subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00370,KO:K08345_1.7.99.4", case=False), "nitrate reductase / nitrite oxidoreductase, alpha subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00370,KO:K08345,KO:K15905_1.7.99.4", case=False), "nitrate reductase / nitrite oxidoreductase, alpha subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00371_1.7.99.4", case=False), " ",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00371,KO:K08346_1.7.99.4", case=False), "nitrate reductase / nitrite oxidoreductase, beta subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00371,KO:K08346,KO:K15906_1.7.99.4", case=False), "nitrate reductase / nitrite oxidoreductase, beta subunit, nitrite oxidase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00372_1.7.99.4", case=False), "assimilatory nitrate reductase catalytic subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00372,KO:K00369_1.7.99.4", case=False), "assimilatory nitrate reductase catalytic subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00372,KO:K02567_1.7.99.4", case=False), "assimilatory nitrate reductase catalytic subunit, periplasmic nitrate reductase NapA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00372,KO:K02567,KO:K00369_1.7.99.4", case=False), "assimilatory nitrate reductase catalytic subunit, periplasmic nitrate reductase NapA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00373", case=False), "nitrate reductase molybdenum cofactor assembly chaperone NarJ/NarW",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00374_1.7.99.4", case=False), "nitrate reductase gamma subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00374,KO:K08347_1.7.99.4", case=False), "nitrate reductase gamma subunit, ",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00376_1.7.2.4", case=False), "nitrous-oxide reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K00376_1.7.99.6", case=False), "nitrous-oxide reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K01425_3.5.1.2", case=False), "glutaminase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K01915_6.3.1.2", case=False), "glutamine synthetase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K01953_6.3.5.4", case=False), "asparagine synthase (glutamine-hydrolysing)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02305", case=False), "nitric oxide reductase subunit C",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02305_1.7.99.7", case=False), "nitric oxide reductase subunit C",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02567_1.7.99.4", case=False), "periplasmic nitrate reductase NapA",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02567,KO:K00372_1.7.99.4", case=False), "periplasmic nitrate reductase NapA, assimilatory nitrate reductase catalytic subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02567,KO:K00372,KO:K00369_1.7.99.4", case=False), "periplasmic nitrate reductase NapA, assimilatory nitrate reductase catalytic subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02568", case=False), "cytochrome c-type protein NapB",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02586_1.18.6.1", case=False), "nitrogenase molybdenum-iron protein alpha chain",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02588_1.18.6.1", case=False), "nitrogenase iron protein NifH",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K02591_1.18.6.1", case=False), "nitrogenase molybdenum-iron protein beta chain",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K03385_1.7.2.2", case=False), "nitrite reductase (cytochrome c-552)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K04013", case=False), "cytochrome c-type protein",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K04014", case=False), "protein",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K04015", case=False), "protein",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K04561_1.7.2.5", case=False), "nitric oxide reductase subunit B",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K05597_3.5.1.38", case=False), "glutamin-(asparagin-)ase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K05601_1.7.99.1", case=False), "hydroxylamine reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K05601,KO:K00368_1.7.99.1,1.7.2.1", case=False), "hydroxylamine reductase, nitrite reductase (NO-forming)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K05601,KO:K15864,KO:K00368_1.7.99.1,1.7.2.1", case=False), "hydroxylamine reductase, nitrite reductase (NO-forming) / hydroxylamine reductase, nitrite reductase (NO-forming)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K08347_1.7.99.4", case=False), "nitrate reductase gamma subunit",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K08361", case=False), "nitrate reductase molybdenum cofactor assembly chaperone NarJ/NarW",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K10534_1.7.1.3", case=False), "nitrate reductase (NAD(P)H)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K10535_1.7.3.4", case=False), "hydroxylamine dehydrogenase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15864", case=False), "nitrite reductase (NO-forming) / hydroxylamine reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15864_1.7.2.1,1.7.99.1", case=False), "nitrite reductase (NO-forming) / hydroxylamine reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15864,KO:K00368_1.7.2.1,1.7.99.1", case=False), "nitrite reductase (NO-forming) / hydroxylamine reductase, nitrite reductase (NO-forming)",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15864,KO:K00368,KO:K05601_1.7.2.1,1.7.99.1", case=False), "nitrite reductase (NO-forming) / hydroxylamine reductase, nitrite reductase (NO-forming), hydroxylamine reductase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15878", case=False), "rieske iron-sulfur protein",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15878,KO:K03890", case=False), "rieske iron-sulfur protein",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15906", case=False), "nitrite oxidase",
                  pd.np.where(Norm_Counts_df.KO_Name.str.contains("KO:K15906,KO:K00371,KO:K08346_1.7.99.4", case=False), "nitrite oxidase, nitrate reductase / nitrite oxidoreductase, beta subunit", "No_name"))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
    #KO_to_gene = pd.read_csv("KO_to_gene.txt", header=0, sep='\t')
    #Norm_Counts_df = pd.merge(Norm_Counts_df, KO_to_gene[["Gene","KO_Name","Description"]], on='KO_Name')
    Norm_Counts_df.to_csv('GeneCountOutput_Foam_Nitro.csv', sep=',',index=False) #write DataFrame to comma separated file (.csv) with file name and FOAM hmm counts
    print("\n"+'\033[93m' + "Printed final file to GeneCountOutput_Foam_Nitro.csv")
    feather.write_dataframe(Norm_Counts_df, "Hmm_hit_counts.feather") #Write to file for reading from R

if __name__ == '__main__':
    main()
