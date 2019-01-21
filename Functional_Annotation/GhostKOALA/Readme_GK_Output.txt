#GhoastKOALA server only takes files under 300MB so file needs to be split to reduce the size.
#Have to double check that a sequence isn't cut off into two files when you do this.
split -b 250M -d protein-sequences.fa

#As of right now there is some blank amino acid sequences. The following warning was produced when running Cogs in anvio.

#WARNING
#===============================================
#63 entries in the sequences table had blank sequences :/ This is related to the
#issue at https://github.com/merenlab/anvio/issues/565. If this is like mid-2018
#and you still get this warning, please find an anvi'o developer and make them
#feel embarrassed. If it is earlier than take this as a simple warning that some
#gene calls in your downstream analyses may have no sequences, and that's OK.
#This is a very minor issue due to on-the-fly addition of Ribosomal RNA gene
#calls to the contigs database, and will likely will not affect anything major.
#This warning will go away when anvi'o can seamlessly work with multiple gene
#callers (which we are looking forward to implement in the future).

#For right now Meren suggests to just remove the blank sequences, with the following script
#https://github.com/merenlab/anvio/issues/943
sed '/^>/ {N; /\n$/d}' protein-sequences2.fa > protein-sequences2_clean.fa

#After running GhoastKOALA you will need to download the files from the webpage and concationate the files 
cat user_ko1.txt user_ko2.txt > user_ko_all.txt
cat user_ko_definition1.txt user_ko_definition2.txt > user_ko_definition_all.txt
cat user.out.top1.gz user.out.top2.gz > user.out.top_all.txt

You can now get count frequences out for each level of Kegg Category with the GK_Count_Frequency.py

#The output of GhostKOALA can be visulaized in iPATH3. Kegg numbers are extracted from GhostKOALA output
awk '{print $2}' user_ko_all.txt | sed '/^[[:space:]]*$/d' > iPalA.txt ; awk 'NR==1{print $1, "color"; next} {print $1, "#ff0000", $2}' iPalA.txt > iPal.txt

#To get GhostKOALAs KOs into Anvio: first we will add the necesary header lines
echo -e "contig\taccession_id" > .temp && cat user_ko_all.txt >> .temp && mv .temp user_ko_all.txt
