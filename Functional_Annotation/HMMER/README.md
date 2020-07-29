After cleaning reads (ran through Kneaddata) reads were merged with flash2 as this might improve hits. 

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
	echo "#!/bin/bash" > ${NEW[$i]}.flash.sh
	echo "##" >> ${NEW[$i]}.flash.sh
	echo "#SBATCH --mem=10G" >> ${NEW[$i]}.flash.sh
	echo "#SBATCH --time=2-18:0:0" >> ${NEW[$i]}.flash.sh
	echo "#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Merged_2018/Error_Out_Files/${NEW[$i]}.Flash.out" >> ${NEW[$i]}.flash.sh
	echo "#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Merged_2018/Error_Out_Files/${NEW[$i]}.Flash.err" >> ${NEW[$i]}.flash.sh
	echo "module load flash2/c41a82e" >> ${NEW[$i]}.flash.sh
	echo "flash2 /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${NEW[$i]}/${NEW[$i]}_L005_R1_001_kneaddata_paired_1.fastq /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${NEW[$i]}/${NEW[$i]}_L005_R1_001_kneaddata_paired_2.fastq -d /share/tearlab/Maga/Jill/CDRF_MetaGenome/Merged_2018/ -o ${NEW[$i]} -r 150 -f 350 -s 35" >> ${NEW[$i]}.flash.sh
	sbatch ${NEW[$i]}.flash.sh
done
```
Statistics on merging can be found in Merged_Stats.txt  

Both merged and forward reads underwent a 6 frame translation with 'sixFrameTranslate.pl'  

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
	echo "#!/bin/bash" > ${NEW[$i]}.6frame.sh
	echo "module load perl-modules/5.18.2" >> ${NEW[$i]}.6frame.sh
	echo "##" >> ${NEW[$i]}.6frame.sh
	echo "#SBATCH --mem=10G" >> ${NEW[$i]}.6frame.sh
	echo "#SBATCH --time=6-18:0:0" >> ${NEW[$i]}.6frame.sh
	echo "#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Protein_2018/Error_Out_Files/${NEW[$i]}.6frame.out" >> ${NEW[$i]}.6frame.sh
	echo "#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Protein_2018/Error_Out_Files/${NEW[$i]}.6frame.err" >> ${NEW[$i]}.6frame.sh
	echo "time perl /share/tearlab/Maga/Jill/bin/sixFrameTranslate.pl -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/Merged_2018/${NEW[$i]}.extendedFrags.fastq -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Protein_2018/Merged_Results/${NEW[$i]}.6frame.faa" >> ${NEW[$i]}.6frame.sh	
	echo "echo "done translating merged reads"" >> ${NEW[$i]}.6frame.sh
	echo "time perl /share/tearlab/Maga/Jill/bin/sixFrameTranslate.pl -i /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${NEW[$i]}/${NEW[$i]}_L005_R1_001_kneaddata_paired_1.fastq -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Protein_2018/Forward_Results/${NEW[$i]}.6frame.faa" >> ${NEW[$i]}.6frame.sh	
	echo "echo "done translating forward reads"" >> ${NEW[$i]}.6frame.sh
	sbatch ${NEW[$i]}.6frame.sh
done
```
hmmscan was used to search the reads against the same subset of nitrogen cycle hmms described in the FOAM_2018 folder.

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
	echo "#!/bin/bash" > NFF.${NEW[$i]}.sh
	echo "##" >> NFF.${NEW[$i]}.sh
	echo "#SBATCH -n 20" >> NFF.${NEW[$i]}.sh
	echo "#SBATCH --mem=40G" >> NFF.${NEW[$i]}.sh
	echo "#SBATCH -p production" >> NFF.${NEW[$i]}.sh
	echo "#SBATCH --time=9-0:0:0" >> NFF.${NEW[$i]}.sh
	echo "#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Foam_2018/Error_Out_Files/${NEW[$i]}.foam.Forward.out" >> NFF.${NEW[$i]}.sh
	echo "#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Foam_2018/Error_Out_Files/${NEW[$i]}.foam.Forward.err" >> NFF.${NEW[$i]}.sh
	echo "module load hmmer/3.1b2" >> NFF.${NEW[$i]}.sh
	echo "cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Foam_2018/" >> NFF.${NEW[$i]}.sh
	echo "#time hmmscan --noali --cpu 30 -E 0.00001 --tblout ./Merged_Results/${NEW[$i]}.Foam.tblout Nitro_Foam.hmm ../Protein_2018/Merged_Results/${NEW[$i]}.6frame.faa" >> NFF.${NEW[$i]}.sh
	echo "#echo "Done Hmmscan on Merged Reads"" >> NFF.${NEW[$i]}.sh
	echo "time hmmscan --noali --cpu 20 -E 0.00001 --tblout ./Forward_Results/${NEW[$i]}.Foam.tblout Nitro_Foam.hmm ../Protein_2018/Forward_Results/${NEW[$i]}.6frame.faa" >> NFF.${NEW[$i]}.sh
	echo "echo "Done Hmmscan on Forward Read"" >> NFF.${NEW[$i]}.sh
	sbatch NFF.${NEW[$i]}.sh
done
```

Microbe Census was run on the reads as well.

```
for f in /share/tearlab/Maga/Jill/CDRF_MetaGenome/Merged_2018/*.extendedFrags.fastq
do
	fname=$(basename $f .extendedFrags.fastq)  
	echo $fname
	echo "#!/bin/bash" > $fname.Mcencus.sh
	echo "##" >> $fname.Mcencus.sh
	echo "#SBATCH -p animal_sciences" >> $fname.Mcencus.sh
	echo "#SBATCH --mem=10G" >> $fname.Mcencus.sh
	echo "#SBATCH --time=2-18:0:0" >> $fname.Mcencus.sh
	echo "#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Nitro_2018/Microbe_Census/Error_Out_Files/${fname}.Microbecensus.out" >> $fname.Mcencus.sh
	echo "#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Nitro_2018/Microbe_Census/Error_Out_Files/${fname}.Microbecensus.err" >> $fname.Mcencus.sh
	echo "run_microbe_census.py -v -t 8 /share/tearlab/Maga/Jill/CDRF_MetaGenome/Merged_2018/${fname}.extendedFrags.fastq /share/tearlab/Maga/Jill/CDRF_MetaGenome/Nitro_2018/Microbe_Census/${fname}_Merged_output.txt" >> $fname.Mcencus.sh
	echo "run_microbe_census.py -v -t 8 /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${fname}/${fname}_L005_R1_001_kneaddata_paired_1.fastq /share/tearlab/Maga/Jill/CDRF_MetaGenome/Nitro_2018/Microbe_Census/${fname}_Forward_output.txt" >> $fname.Mcencus.sh
sbatch $fname.Mcencus.sh
done
```
The results for forward reads are found in FwdReads_MicrobeCensus_Results.txt 

To get counts of these genes that can be read into R run the follow:

```
 python FOAM_hmmscan_to_GeneCount.py -i ./ -r Merged_Results.txt --recA GeneCountOutput_Nitro.csv -s Merge_Stats.txt
```
This produces `GeneCountOutput_Foam_Nitro.csv` that contains the information for graphing. 
