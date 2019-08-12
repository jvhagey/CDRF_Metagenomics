After running Prokka:  
`prokka --cpus 20 --metagenome --outdir ../Prokka/Output_redo/ --prefix Standard fixed.contigsV5_1000.fa`

We will take a look at the GFF file made by Prokka:  
`less -S ./Output/Standard.gff`
You will notice that the first lines in the GFF file show the annotated sequence regions.  
To skip these and get directly to the annotations you can do:  
`grep -v "^#" ./Output/Standard.gff | less -S`  

Some genes in the dataset should now contain annotations from several databases, such as enzyme comission and COG (Clusters of Orthologous Groups) identifiers. In the downstream analyses we will quantify and compare the abundance of enzymes and metabolic pathways, as well as COGs in the different samples. To do this, we will first extract lists of the genes with enzyme and COG IDs from the GFF file that was produced by PROKKA.

First, we find lines containing enzyme numbers using pattern matching with grep and then reformat the output using a combination of cut and sed:
`grep "eC_number=" ./Output/Standard.gff | cut -f9 | cut -f1,2 -d ’;’| sed ’s/ID=//g’| sed ’s/;eC_number=/\t/g’ > Standard.ec`
Then we extract COG identifiers:
`egrep "COG[0-9]{4}" ./Output/Standard.gff | cut -f9 | sed 's/.\+COG\([0-9]\+\).\+;locus_tag=\(GANJLKBE_[0-9]\+\);.\+/\2\tCOG\1/g' > Standard.cog`

We will now have two tab deliminated files that have gene calls followed by the E.C. number or COGID.
Next we downloaded the most recent data (2014 update) from NCBI's COG web page https://www.ncbi.nlm.nih.gov/COG/
The Standard.ec and Standard.cog file can be parced into something that R can used with XXX.py and Prokka_Cog_Counts.py respectively. 

JGI also has a list of Cogs: https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=FindFunctions&page=cogid2cat

MinPath
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/MinPath/
python2 /share/magalab/bin/MinPath/MinPath1.4.py -ko ../GhostKOALA/user_ko_all.txt -report Minpath_ko.report -details Minpath_ko.details
