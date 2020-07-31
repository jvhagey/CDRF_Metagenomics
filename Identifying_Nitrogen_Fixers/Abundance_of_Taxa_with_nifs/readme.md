## For the nifH genes found in the same contig with nifD and nifK

### Identifying Taxa with nifHDK in the same contig
Run `python run_all_works2.py -i ./ -g "nifH|nifK|nifD|COG2710" -t family`. The output of this is found in the file `Families_with_Nif_Genes.txt`. For now we will look just for taxa with nifH `python run_all_works2.py -i ./ -g "nifH" -t family -o some_name_you_pick.txt` and you can swap family out for other levels. The output of this is found in `Taxa_with_nifH.txt`. Now we will move on to getting the gene calls for nifH and placing them on a tree. See the nifH_Tree folder. 


## For the nifH genes found in a contig with no other nif genes

### Identifying Taxa with nif in separate contigs

### Testing Differential Abundance

I first extracted the families of interest (Families_with_Nif_Genes.txt) from the output of `kaiju-addTaxonNames` files.

```
for f in /share/tearlab/Maga/Jill/CDRF_MetaGenome/Kaiju_2019/*_kaiju-names_out.tsv
do
	fname=$(basename $f _kaiju-names_out.tsv)
	#grep to find families of interest in each of the kaiju folders
	time grep -f Families_of_interest.txt ${fname}_kaiju-names_out.tsv >> ${fname}_kaiju-names_out_subset.tsv
done
```
Note: that in the script above `Families_of_interest.txt` is a file with the name on a family on each line. 
Script `Kaiju_to_otutab.py` was written to parse the output of `kaiju-addTaxonNames` into a tax_tab.txt and otu_tab.txt that can be read into R to create a phyloseq object. For some reason I haven't sorted out yet it isn't running at as an executable, but will run by line. Shoot me an email if you want some explaination or it fixed for your own purposes.

```
time python Kaiju_to_otutab.py -f _kaiju-names_out_subset.tsv
```
