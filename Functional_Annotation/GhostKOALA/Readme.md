GhostKOALA  
First, run extract gene calls from Anvi'o database.
```
#!/bin/bash
##
#SBATCH --mem=20G
#SBATCH --time=3-18:0:0
#SBATCH -n 3
#SBATCH -p production
#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/Error_Out_Files/Get_AA.out
#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/Error_Out_Files/Get_AA.err

source activate Anvio5.5

#get gene calls to run GhostKOALA
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/
time anvi-get-sequences-for-gene-calls --get-aa-sequences -c fixed.contigsV5.5.db -o ./GhostKOALA/protein-sequences.fa
echo "Done getting AA calls"
#adding genecall_ to the end of reach gene call so that the GhoastKOALA server isn't mad at us
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/GhostKOALA/
sed -i 's/>/>genecall_/g' protein-sequences.fa
echo "Done Changing gene call names"
```

GhostKOALA server only takes files under 300MB so file needs to be split to reduce the size. Note, you will have to double check that a sequence isn't cut off into two files when you do this.  
`split -b 250M -d protein-sequences.fa`  

This gives two files x01 and x00, they were renamed to protein-sequences_1.fa and protein-sequences_2.fa.  

After, running GhoastKOALA searching the KEGG database: c_family_euk+genus_prok+viruses you will need to download the files from the webpage and concationate the files. The files are then concatenated back into one file.  

```
cat user_ko_1.txt user_ko_2.txt > user_ko_all.txt
cat user_ko_definition_1.txt user_ko_definition_2.txt > user_ko_definition_all.txt
cat user.out.top_1.gz user.out.top_2.gz > user.out.top_all.txt
cat kegg_taxonomy_1.list kegg_taxonomy_2.list > kegg_taxonomy_all.list
```

You can now get count frequences out for each level of Kegg Category with the GK_Count_Frequency.py

The output of GhostKOALA can be visulaized in iPATH3. Kegg numbers are extracted from GhostKOALA output.  
`awk '{print $2}' user_ko_all.txt | sed '/^[[:space:]]*$/d' | sed "1iKO" > iPalA.txt; awk 'NR==1{print $1, "color", "width"; next} {print $1, "#ff0000", "W10" $2}' iPalA.txt > iPal.txt`

To get GhostKOALAs KOs into Anvio you can follow [this](http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/) tutorial. First, we will add the necesary header lines.   
`echo -e "contig\taccession_id" > .temp && cat user_ko_all.txt >> .temp && mv .temp user_ko_all.txt`

Download the parser for GhostKOALA  
`git clone https://github.com/edgraham/GhostKoalaParser.git`  

Now we can import GhostKOALA annotations into Anvio database.  

```
time KEGG-to-anvio --KeggDB /share/magalab/bin/GhostKOALA_to_Anvi/GhostKoalaParser/samples/KO_Orthology_ko00001.txt -i user_ko_all.txt -o KeggAnnotations-AnviImportable.txt
#importing
time anvi-import-functions -c ../C5_V5.5/fixed.contigsV5.5.db -i KeggAnnotations-AnviImportable.txt
```
