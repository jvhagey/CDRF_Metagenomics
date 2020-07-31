# Taxonamy Assignment with Kaiju  
Building database for Kaiju first  

```
#!/bin/bash
##
#SBATCH -p gc128,gc512
#SBATCH --mem=450G
#SBATCH --time=4-18:0:0
#SBATCH -o /share/tearlab/Maga/Jill/Kaijudb_nr2.out
#SBATCH -e /share/tearlab/Maga/Jill/Kaijudb_nr2.err

#first open folder where I want to the db to be created. 
#Non-redundant protein database Archaea, Bacteria, Fungi, Microbial Eukaryotes and Viruses.
cd /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/
time /share/tearlab/Maga/Jill/bin/kaiju/bin/makeDB.sh -e

cd /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/
time /share/tearlab/Maga/Jill/bin/kaiju/bin/mkfmi kaiju_db_nr_euk
```

Running Kaiju  

```
#!/bin/bash

NEW[0]=1-10_S26
NEW[1]=1-12_S27
NEW[2]=1-1_S25
NEW[3]=5-11_S29
NEW[4]=5-14_S30
NEW[5]=5-6_S28
NEW[6]=6-1_S31
NEW[7]=6-5_S32
NEW[8]=6-7_S33
NEW[9]=8-10_S36
NEW[10]=8-8_S34
NEW[11]=8-9_S35

for i in "${!NEW[@]}"
do
	echo ${NEW[$i]}
	echo "#!/bin/bash" > ${NEW[$i]}.Kaiju.sh
	echo "##" >> ${NEW[$i]}.Kaiju.sh
	echo "#SBATCH -p gc64,gc128,gc512,gc256" >> ${NEW[$i]}.Kaiju.sh
	echo "#SBATCH --mem=120G" >> ${NEW[$i]}.Kaiju.sh
	echo "#SBATCH --time=6-18:0:0" >> ${NEW[$i]}.Kaiju.sh
	echo "#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/Error_Out_Files/${NEW[$i]}.Kaiju.out" >> ${NEW[$i]}.Kaiju.sh
	echo "#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/Error_Out_Files/${NEW[$i]}.Kaiju.err" >> ${NEW[$i]}.Kaiju.sh
	echo "#Note this script (line 27) ran out of memory with 50G and 6 threads. Needed 120G and 10 threads" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/kaiju -z 10 -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -f /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/kaiju_db_nr_euk.fmi -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${NEW[$i]}/${NEW[$i]}_L005_R1_001_kneaddata_paired_1.fastq -j /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${NEW[$i]}/${NEW[$i]}_L005_R1_001_kneaddata_paired_2.fastq -o ${NEW[$i]}_kaiju_out" >> ${NEW[$i]}.Kaiju.sh
	echo "#will make report down to genus level of all classified taxa in report" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/kaijuReport -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -n /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/names.dmp -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/${NEW[$i]}_kaiju_out -r genus -o ${NEW[$i]}_kaiju_out_genus.summary" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/kaijuReport -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -n /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/names.dmp -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/${NEW[$i]}_kaiju_out -r species -o ${NEW[$i]}_kaiju_out_species.summary" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/kaijuReport -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -n /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/names.dmp -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/${NEW[$i]}_kaiju_out -r class -o ${NEW[$i]}_kaiju_out_class.summary" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/kaijuReport -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -n /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/names.dmp -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/${NEW[$i]}_kaiju_out -r order -o ${NEW[$i]}_kaiju_out_order.summary" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/kaijuReport -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -n /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/names.dmp -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/${NEW[$i]}_kaiju_out -r phylum -o ${NEW[$i]}_kaiju_out_phylum.summary" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/kaijuReport -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -n /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/names.dmp -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/${NEW[$i]}_kaiju_out -r family -o ${NEW[$i]}_kaiju_out_family.summary" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/addTaxonNames -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -n /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/names.dmp -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/${NEW[$i]}_kaiju_out -p -o ${NEW[$i]}_kaiju-names_out" >> ${NEW[$i]}.Kaiju.sh
	echo "time /share/tearlab/Maga/Jill/bin/kaiju/bin/kaiju2krona -t /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/nodes.dmp -n /share/tearlab/Maga/Jill/bin/kaiju/kaijudb_nr/names.dmp -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2018/Kaiju_Output/${NEW[$i]}_kaiju-names_out -o ${NEW[$i]}_kaiju-names_out.krona" >> ${NEW[$i]}.Kaiju.sh
	sbatch ${NEW[$i]}.Kaiju.sh
done
```
Script `Kaiju_to_R.py` was written to parse kaiju output .summary files into a .csv file that can be read into R. With Kaiju's newest update `Kaiju_v1.7.2_to_R.py` was created to handled Kaiju's new output, which was significantly cleaner. Both have the following options:

```
usage: Kaiju_v1.7.2_to_R.py [-h] -t TAXLEVEL -i INDIRNAME

Parses Kaiju .summary files into an R readable .csv file

optional arguments:
  -h, --help            show this help message and exit
  -t TAXLEVEL, --taxonomic-level TAXLEVEL
                        Pick One: phylum, class, order, family, genus, species
                        (required)
  -i INDIRNAME, --indir INDIRNAME
                        Input directory name were .summary files are found
                        (required)
```
