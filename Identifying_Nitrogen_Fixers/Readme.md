Once the Anvi'o database was constructed tables of functions, taxonomy and hmm hits were exported. 
```
#!/bin/bash
##
#SBATCH -p production
#SBATCH --mem=10G
#SBATCH -n 2
#SBATCH --time=10:0:0
#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/Error_Out_Files/Anvi_export.out
#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/Error_Out_Files/Anvi_export.err

aklog
source activate Anvio5.5

cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/

#checking list of functions available
anvi-export-functions -c fixed.contigsV5.5.db -l

time anvi-export-functions -c fixed.contigsV5.5.db -o ./Anvi_Tables/Functions.txt
echo "Done exporting Functions"

#Get list of tables available
anvi-export-table -l fixed.contigsV5.5.db

time anvi-export-table --table genes_taxonomy -o ./Anvi_Tables/genes_taxonomy.txt fixed.contigsV5.5.db
echo "Done genes_taxonomy"
time anvi-export-table --table hmm_hits -o ./Anvi_Tables/hmm_hits.txt fixed.contigsV5.5.db
echo "Done hmm_hits"
time anvi-export-table --table hmm_hits_info -o ./Anvi_Tables/hmm_hits_info.txt fixed.contigsV5.5.db
echo "Done hmm_hits_info"
time anvi-export-table --table gene_functions -o ./Anvi_Tables/gene_functions.txt fixed.contigsV5.5.db
echo "Done gene_functions"
time anvi-export-table --table taxon_names -o ./Anvi_Tables/taxon_names.txt fixed.contigsV5.5.db
echo "Done taxon_names"
time anvi-export-table --table genes_in_contigs -o ./Anvi_Tables/genes_in_contigs.txt fixed.contigsV5.5.db
echo "Done genes_in_contigs"

time anvi-export-splits_taxonomy -c fixed.contigsV5.5.db -o ./Anvi_Tables/splits_taxonomy.txt

time anvi-export-table --table splits_taxonomy -o ./Anvi_Tables/splits_taxonomy.txt fixed.contigsV5.5.db
echo "Done splits_taxonomy"
time anvi-export-table --table hmm_hits_in_splits -o ./Anvi_Tables/hmm_hits_in_splits.txt fixed.contigsV5.5.db
echo "Done hmm_hits_in_splits"
time anvi-export-table --table genes_in_splits -o ./Anvi_Tables/genes_in_splits.txt fixed.contigsV5.5.db
echo "Done genes_in_splits"
```

Next we parsed the Anvi'o database tables to find contigs of interest that have hits against genes in the nitrogen cycle. First, we identified contigs that had hits with nifHDK in them and were identified this way with Foam_HMMs, GhostKOALA and COGs. Next, we identified contigs that had only one of these nif genes. Bins will be searched for contigs containing nifHDK as well as those having individual contigs for each gene.

```
python ../run_all_works2.py -i ./ --get-HDK-contigs -o Allmethods_HDK_contigs.txt
python ../run_all_works2.py -i ./ -g "nifH" -o Allmethods_nifH_contigs.txt
python ../run_all_works2.py -i ./ -g "nifK|COG2710|nifD" -o Allmethods_nifKD_contigs.txt
```

Once we have a list of contigs of interest we binned contigs with MetaBat2.

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

After running MetaBAT2, the list of contigs of interest were searched for in the bins. And these bins were placed in a new folder.  
The first step of this is taking the list of contigs names in a text file created by using the parser `python -i ./ `.  Note that I used seqtk, but you can also use qiime1 script `filter_fasta.py` see [here](http://qiime.org/scripts/filter_fasta.html).  

```
#!/bin/bash
##
#SBATCH --mem=5G
#SBATCH -p production
#SBATCH -n 3
#SBATCH --time=6-18:0:0
#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/Error_Out_Files/metabat_anvio_contigs_nif.out
#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/Error_Out_Files/metabat_anvio_contigs_nif.err

module load seqtk/1.3

#using seqk search fasta file for specific contigs.
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/
time seqtk subseq ../C5_V5.1/fixed.contigsV5_1000.fa ../SourMash/GK_NifD_contigs.txt > ../SourMash_2019/contigs_GK_NifD.fa

```
Then pass grep a file that has contigs in it and search all bins in a folder to get list of file name (AKA the bin) with the contigs of interest in.  
```
grep -H -f contigs_GK_NifH.fa ../../MetaBat/C5_Anvio/*.fa > bins_with_GKnifH_Contigs.txt
```
Next, parse the list of bins to just get the unique bins that contained these contigs.  

```
sed 's/:>c_[0-9]\+$//' bins_with_GKnifH_Contigs.txt | sed 's/.*C5_Anvio//' | sed 's/^/C5_Anvio/'| uniq > bins_with_GKnifH_Contigs2.txt
```
Lastly, you can make a new directory and copy the bins you identified into a new folder and then run Sourmash on that folder of bins of interest to further confirm species of interest that might be able to nitrogen fix.  

```
#Make directorys
mkdir MetaBat_bins_NifH
mkdir MetaBat_bins_NifD
mkdir MetaBat_bins_NifK

#copy the fasta files into new folder.
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/MetaBat/C5_Anvio/

cat ../../Anvio/SourMash/bins_with_GKnifH_Contigs2.txt | xargs cp -t ../../Anvio/SourMash/MetaBat_bins_NifH/


```



Sourmash was run on each bin to assign taxonomy. 

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
[Clostridium sp. Uncultured A]() Accession BioProject #PRJEB10915 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=FMGR01)  
[Ruminococcaceae bacterium YAD300](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB16616) Assembly Accession #GCA_900107185.1 BioProject #PRJEB16616 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/FNPA01?display=contigs)  
[Sarcina sp. DSM 11001](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB16054) Accession #FNFZ00000000 BioProject #PRJEB16054 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=FNFZ01)  
[Clostridium sp. Uncultured B]() Accession #FMGT01000000 BioProject #PRJEB10915 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=FMGT01)  
[Eubacterium sp. AB3007](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA217194) Accession #JIAD00000000 BioProject #PRJNA217194 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=JIAD01)  
[Erysipelotrichaceae bacterium NK3D112](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA223066) Accession #JNJQ00000000 BioProject #PRJNA223066 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=JNJQ01)  
[Clostridium celatum DSM 1785](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA30375) Accession AMEZ00000000 BioProject #PRJNA30375 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=AMEZ01)  
[Anaerostipes sp. 992a](https://www.ncbi.nlm.nih.gov/nuccore/MJIG00000000) Accession #MJIG00000000 BioProject #PRJNA341691 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=MJIG01)  
[Prevotella bryantii](https://www.ncbi.nlm.nih.gov/nuccore/NPJF00000000) Accession #NPJF00000000 BioProject #PRJNA396820 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=NPJF01)  
[Lachnobacterium bovis](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB16629) Accession #FNPG00000000 BioProject #PRJEB16629 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=FNPG01)  
[Lachnobacterium bacterium](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB19693) Accession #FUZG00000000 BioProject #PRJEB19693 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=FUZG01)  
[Olsenella umbonata](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB15778) Accession #FNWT00000000 BioProject #PRJEB15778 [Sequence](https://www.ncbi.nlm.nih.gov/Traces/wgs/CYZX01?val=FNWT01)  
