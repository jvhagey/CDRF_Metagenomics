## For the nifH genes found in the same contig with nifD and nifK

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
time anvi-get-sequences-for-gene-calls -c fixed.contigsV5.5.db --get-aa-sequences --gene-caller-ids 1299485,1401422,264012 -o nifH_AA_genecalls.fasta
``` 
Now that we have sequences we will align them to a reference dataset of prealigned nifH sequences from the [the Zehr Lab](https://www.jzehrlab.com/nifh) that are in the file `CART_Test_Atlantic.fasta`. You will need to follow the CART button to get to this file. Not that to my gene call file, I added a couple sequences as positive controls to make sure the script worked correctly.

```
source activate MAFFT

mafft --add nifH_AA_genecalls.fasta --reorder CART_Test_Atlantic.fasta > Atlantic_withMyNifSeqs.fasta
```

`NifH_Clusters.py` is a script to determine where putative nitrogen fixing sequences sit on a NifH tree. This script came from [the Zehr Lab](https://www.jzehrlab.com/nifh). 
The script is descripted in [Frank et al. 2016](https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/1758-2229.12455)

```
python NifH_Clusters.py Atlantic_withMyNifSeqs 45
```

To get the taxonomy of these gene calls run the following. 

```

python anvi-get-taxonomy-for-genecall.py -i ./ -g "139447" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "264012" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "1299485" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "1401422" -t species

```
I added these manually to the fasta file.

This assigned the sequences to the following clusters/subclusters:

```
>139447 | Lachnospiraceae bacterium XPB1003 main cluster = 3 subcluster = 3I
>264012 | Prevotella bryantii main cluster = 3 subcluster = 3I
>1299485 | Unknown_Methanobrevibacter main cluster = 4 subcluster = 4B
>1401422 | Unknown_Methanobrevibacter main cluster = 4 subcluster = 4B
```
## For the nifH genes found in a contig with no other nif genes, but are in a bin with nifDK

Again we will get the gene calls for nifH that are in a contig that is in a bin with another contig(s) that have nifDK.

```
python get_gene_calls_from_contig.py -i ./ -c "c_000000171271|c_000002173983|c_000002488619|c_000002799447|c_000000491992|c_000000529426|c_000001311790|c_000001523212|c_000003416518|c_000005865692" -tf "TIGR01287" -o nifH_genecalls_solocontig_df.txt
```

Next we get the actual fasta sequence for these gene calls.

```
time anvi-get-sequences-for-gene-calls -c fixed.contigsV5.5.db --get-aa-sequences --gene-caller-ids 182510,472009,1154496,64175,760920,196015,196014,543781,962682,760926,863085 -o nifH_solo_AA_genecalls.fasta
``` 

Next we grab the taxonomy of these gene calls by running the following. 

```
python anvi-get-taxonomy-for-genecall.py -i ./ -g "182510" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "472009" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "1154496" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "64175" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "760920" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "196015" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "196014" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "543781" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "962682" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "760926" -t species
python anvi-get-taxonomy-for-genecall.py -i ./ -g "863085" -t species
```

I added these manually to the fasta file then aligned these sequneces to nifH clusters and subclusters reference sequences.

```
source activate MAFFT

mafft --add nifH_solo_AA_genecalls.fasta --reorder CART_Test_Atlantic.fasta > Atlantic_withMyNifSoloSeqs.fasta
```

This assigned the sequences to the following clusters/subclusters:

```
python NifH_Clusters.py Atlantic_withMyNifSoloSeqs 45
```

```
>182510 | Lachnospiraceae bacterium XPB1003 main cluster = 3 subcluster = 3I

```
