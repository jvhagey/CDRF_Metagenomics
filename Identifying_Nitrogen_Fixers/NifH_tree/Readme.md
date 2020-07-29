

The fasta file for gene calls for each gene in the text file were extracted from the file of all gene calls from the anivo database. 

```
#!/bin/bash
##
#SBATCH --mem=5G
#SBATCH -p production
#SBATCH -n 3
#SBATCH --time=6-18:0:0
#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/Error_Out_Files/metabat_anvio_contigs_nif.out
#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/Error_Out_Files/metabat_anvio_contigs_nif.err

aklog
module load seqtk/1.3
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/Gene_calls/

for f in /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/Gene_calls/*_genecalls.txt
do
	fname=$(basename $f _genecalls.txt)
	#using seqk search fasta file for specific contigs.
	time seqtk subseq ../C5_V5.5_1000_gene-calls.fa ${fname}_genecalls.txt > ${fname}_genecalls.fa
done
```



`NifH_Clusters.py` is a script to determine where putative nitrogen fixing sequences sit on a NifH tree. This script came from [the Zehr Lab](https://www.jzehrlab.com/nifh). 
The script is descripted in [Frank et al. 2016](https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/1758-2229.12455)
