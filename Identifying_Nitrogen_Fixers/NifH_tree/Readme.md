Here we want to get the gene call Ids to be able to pull the fasta files for each genne call from the Anvio'd database. To get gene call IDs is a bit of an "easter egg" run the following command `python run_all_works2.py -i ./ -g "nifK|COG2710|nifD" -t family -o test.txt`. You will be prompted with questions and you will need to answer `contigs`, `no` and `per gene`. This part of the script isn't set up to get a particular set of genes yet :(. This will create a bunch of files that have gene calls for each gene. For example, the nifH file will have gene calls IDs for genes that were identified as nifH by COGs, TIGRFAM HMMs, FOAM HMMs and GhoastKOALA. 

Next we will get a fasta file of all the gene calls in the Anvi'o database.

```
anvi-get-sequences-for-gene-calls -c fixed.contigsV5.5.db --get-aa-sequences -o C5_V5.5_1000_gene-calls_AA.fa
anvi-get-sequences-for-gene-calls -c fixed.contigsV5.5.db -o C5_V5.5_1000_gene-calls.fa

```

Now that we have a text file with gene call IDs we can search the fasta file of gene calls from Anvi'o to extract the sequences we want. 

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

You can also run something like this.

```
time anvi-get-sequences-for-gene-calls -c fixed.contigsV5.5.db --gene-caller-ids 1299485,1401422,264012 -o nifH_DNA_genecalls.fasta
``` 
Now that we have sequences we will align them to a reference dataset of prealigned nifH sequences from the [the Zehr Lab](https://www.jzehrlab.com/nifh) that are in the file `CART_Test_Atlantic.fasta`. You will need to follow the CART button to get to this file. Not that to my gene call file, I added a couple sequences as positive controls to make sure the script worked correctly.

```
module load mafft/7.453

mafft --add nifH_AA_genecalls_notaligned.fasta --reorder CART_Test_Atlantic.fasta > Atlantic_withMyNifSeqs.fasta
```

`NifH_Clusters.py` is a script to determine where putative nitrogen fixing sequences sit on a NifH tree. This script came from [the Zehr Lab](https://www.jzehrlab.com/nifh). 
The script is descripted in [Frank et al. 2016](https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/1758-2229.12455)

```
python NifH_Clusters.py Atlantic_withMyNifSeqs.fasta 45
```
