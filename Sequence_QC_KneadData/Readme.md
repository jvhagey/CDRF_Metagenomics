#Quality and Trimming

The following script was used to analyze metagenomic sequencing of CDRF fecal microbiota  

4 farms with 4 samples from each farms  
Raw Files found at /share/tearlab/Maga/Jill/CDRF_MetaGenome/Raw/  
Each file was trimmed and host & PhiX174 reads removed to ensure quality  
```
for f in /share/tearlab/Maga/Jill/CDRF_MetaGenome/Raw/*_L005_R1_001.fastq.gz
do
	fname=$(basename $f _L005_R1_001.fastq.gz) 
	echo $fname
	echo "#!/bin/bash" > $fname.kneaddata.sh
	echo "##" >> $fname.kneaddata.sh
	echo "#SBATCH -p gc64" >> $fname.kneaddata.sh
	echo "#SBATCH --mem=60G" >> $fname.kneaddata.sh
	echo "#SBATCH --time=4-18:0:0" >> $fname.kneaddata.sh
	echo "#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${fname}.kneaddata.out" >> $fname.kneaddata.sh
	echo "#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${fname}.kneaddata.err" >> $fname.kneaddata.sh
	echo "module load bowtie2/2.2.8" >> $fname.kneaddata.sh
	echo "module load trimmomatic/0.33" >> $fname.kneaddata.sh
	echo "module load fastqc/v0.11.5" >> $fname.kneaddata.sh
	echo "/home/jvhagey/.local/bin/kneaddata --input /share/tearlab/Maga/Jill/CDRF_MetaGenome/Raw/${fname}_L005_R1_001.fastq.gz --input /share/tearlab/Maga/Jill/CDRF_MetaGenome/Raw/${fname}_L005_R2_001.fastq.gz -t 5 -db /share/tearlab/Maga/Jill/bin/Fastq_Screen_Index/Cow_Genome_Index/ -db /share/tearlab/Maga/Jill/bin/Fastq_Screen_Index/Phage_Index/ --output /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/${fname}/ --trimmomatic-options="ILLUMINACLIP:/share/tearlab/Maga/Jill/adapter.fasta:2:30:10" --trimmomatic-options="HEADCROP:10" --trimmomatic-options="LEADING:20" --trimmomatic-options="TRAILING:20" --trimmomatic-options="SLIDINGWINDOW:4:20" --trimmomatic-options="MINLEN:100" --run-fastqc-start --run-fastqc-end" >> $fname.kneaddata.sh
	sbatch $fname.kneaddata.sh
done
```
Note First need to make the indexes from NCBI sequences

```
#!/bin/bash
##
#SBATCH -p gc64
#SBATCH --mem=60G
#SBATCH --time=4-18:0:0
#SBATCH -o /share/tearlab/Maga/Jill/Kneaddata_DB.out
#SBATCH -e /share/tearlab/Maga/Jill/Kneaddata_DB.err
module load bowtie2/2.2.8

#downloaded from https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=DAAA02#contigs
#unzip file 
#gunzip /share/tearlab/Maga/Jill/bin/Fastq_Screen_Index/Cow_Genome/Cow_Genome.fsa.gz
#making bowtie index
time bowtie2-build -f /share/tearlab/Maga/Jill/bin/Fastq_Screen_Index/Cow_Genome/Cow_Genome.fsa /share/tearlab/Maga/Jill/bin/Fastq_Screen_Index/Cow_Genome/Cow_Whole_Genome
#Download from http://www.ncbi.nlm.nih.gov/nuccore/9626372?report=fasta
time bowtie2-build -f /share/tearlab/Maga/Jill/16s/Fastq_screen/phage_phiX174_genome.fasta Phage
```
To get trimming stats:

```
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/KneadData/Error_Out_Files/
grep "Initial number of reads" *.kneaddata.out
grep -A 1 "Running Trimmomatic" *.kneaddata.out
grep "_paired_1.fastq" *.kneaddata.out
```

