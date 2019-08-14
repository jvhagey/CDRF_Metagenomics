Once the Anvi'o database was constructed and parsed to find contigs of interest contigs were binned with MetaBat2

Prior to running Metabat you need to calculate the depth of coverage

```
#Running MetaBAT2 on the "fixed" contigs that were formatted for use with Anvi'o
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/
time jgi_summarize_bam_contig_depths --outputDepth ./MetaBat/depth.txt ./Anvio/C5_Mapping/*.sorted.bam
echo "done calculating depth"
```
Next, MetaBAT2 is run. 
```
#!/bin/bash
##
#SBATCH --mem=500G
#SBATCH -p assembly
#SBATCH -n 25
#SBATCH --time=6-18:0:0
#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/MetaBat/Error_Out_Files/metabat_binning_more.out
#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/MetaBat/Error_Out_Files/metabat_binning_more.err

#needed 500G mem
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5/

time /software/metabat/2.12.1/static/metabat2 -t 25 -i fixed.contigs.fa -a ../../MetaBat/depth.txt -o ../../MetaBat/C5_bin_1500/C5_bin_more --unbinned -v --saveCls --seed 125 --minContig 1500
echo "done Metabat 1500"
```

After running MetaBAT2, the list of contigs of interest were searched for in the bins. And these bins were placed in a new folder. Sourmash was thing run on each bin to assign taxonomy. 

```
#!/bin/bash
##
#SBATCH -p production
#SBATCH -n 1
#SBATCH --time=10:0:0
#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/Error_Out_Files/gather_loop_metabat.out
#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/Error_Out_Files/gather_loop_metabat.err

aklog
source activate sourmash2

cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/

for f in /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/MetaBat_bins/*.fa
do
	fname=$(basename $f .fa )
	time sourmash compute --scaled 2000 ./MetaBat_bins/${fname}.fa -o ./Sig_Files/${fname}.sig -k 21,31,51 >> ${fname}.out
	time sourmash gather ./Sig_Files/${fname}.sig -k 21 /share/magalab/bin/SourMash_DB/Genbank/genbank-d2-k21.sbt.json -o ${fname}.csv --save-matches ./Matches/${fname}_matchsigs.csv >> ${fname}.out
done
```

Once accession numbers had been identified the full genome were downloaded.  
[Clostridium saudiense](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB4942) Accession #CBYM010000000 BioProject #PRJEB4942 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CBYM01?display=contigs)  
[Turicibacter sanguinis PC909](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA42765) Accession #ADMN00000000 BioProject #PRJNA42765 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/ADMN01?display=contigs)  
[Romboutsia timonensis](https://www.ncbi.nlm.nih.gov/bioproject/324389) Accession #FNMT00000000 BioProject #PRJEB14233 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/FNMT01?display=contigs)  
[Clostridium disporicum](https://www.ncbi.nlm.nih.gov/assembly/GCF_001405015.1/#/def) Assembly Accession GCF_001405015.1 BioProject #PRJEB10915 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?display=contigs)  
[Methanobrevibacter sp. YE315](https://www.ncbi.nlm.nih.gov/assembly/GCA_001548675.1) Assembly Accession GCA_001548675.1 BioProject #PRJNA273773 [Sequence](https://www.ncbi.nlm.nih.gov/nuccore/CP010834.1)  
[Methanobrevibacter millerae strain SM9](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA49589) Accession BioProject #PRJNA49589 [Sequence](https://www.ncbi.nlm.nih.gov/nuccore/CP011266)  
[Clostridium sp. Uncultured A]()[Sequence]()
[Ruminococcaceae bacterium YAD300](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB16616) Assembly Accession #GCA_900107185.1 BioProject #PRJEB16616 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/FNPA01?display=contigs)  


