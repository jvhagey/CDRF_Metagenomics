Currently, the arguments accepted for `Anvi_table_parser.py` are:

```
usage: Anvi_table_parser.py [-h] [-t TAXALEVEL] [-l FOAMLEVEL]
                  Identify_Nif_Taxonom              [--print-contigs] [--find-nifHDK-in-contigs]
                                [-st SPECIFICTAXA] -i INDIRNAME [-o OUTFILE]
                                [-g GENE] -s SOURCE

Parses anvi-export-table --table files

optional arguments:
  -h, --help            show this help message and exit
  -t TAXALEVEL, --taxa-level TAXALEVEL
                        Pick taxa level you want reported
  -l FOAMLEVEL, --foam-level FOAMLEVEL
                        For list of foam level hits use '-l list' as arguement
                        default none.
  --print-contigs       Print contigs with your choose combination of source,
                        taxa and gene.
  --find-nifHDK-in-contigs
                        Print contigs with nifHDKENB genes.
  -st SPECIFICTAXA, --specific-taxa SPECIFICTAXA
                        Give name of specific taxa you are interested in. You
                        need to pass -t for the appropriate level as well.
  -i INDIRNAME, --indir INDIRNAME
                        Input directory name were taxon_names, genes_taxonomy
                        and hmm_hits .txt files are found (required)
  -o OUTFILE, -outfile OUTFILE
                        Name of file to write output too.
  -g GENE, --gene GENE  Name of gene to search in dataframe for. You can pass
                        multiple genes by separating them by a pipe. Ex:
                        'geneA|geneB'
  -s SOURCE, --source SOURCE
                        Indicate source, for list of sources use '-s list' as
                        arguement default none.
 ```
 
 To get a list of sources use `python parser.py -i ./ -s "list"`. The sources for this Anvi'o database are:
 
 ```
 ['Campbell_et_al', 'Rinke_et_al', 'Foam_Nitro', 'Nif_Hmms', 'Ribosomal_RNAs', 'KeggGhostKoala', 'COG']
```

If you want a full dataframe for a particular source use `python parser.py -i ./ -s "Foam_Nitro"`. This will output the full dataframe
and will prompt you to see if you want to save it or not. If you want only the dataframe for a particular gene use 
`python parser.py -i ./ -s "Foam_Nitro" -g "niH"`. If you add the taxonomy flag like `python parser.py -i ./ -s "Foam_Nitro" -g "niH" -t family`
to this you will also be prompted to see if you want counts of how many times a taxonomy level is seen with that gene or if you just want a 
list of unique taxa that have that gene. Using the taxonomy flag without a gene will return the entire dataframe and after answering if you want 
the file saved you will be asked how you want to view the taxonomic data. If you plan on printing out a dataframe or list you will have to 
supply a file name. If you add `--print-contigs` in combination with any of the above arguments will give you a text file of these
contigs, just dont forget to set a file name with the `-o` argument. Using the `--find-nifHDK-in-contigs` flag with any of the sources and 
a taxonomy will tell you what taxa, at the level you choose, have contigs that have hits with those genes. 

If you want to find taxa or contigs identified to have nifHDK genes or a specific gene that were found by multiple methods use the script `Find_common_contigs.py`

Currently, the arguments accepted for `Find_common_contigs.py` are:

```
usage: Find_common_contigs.py [-h] -i INDIRNAME [-t TAXALEVEL] [--print-taxa]
                              [--get-taxa-contigs] [-g GENE]
                              [--get-gene-contigs]

Parses anvi-export-table --table files

optional arguments:
  -h, --help            show this help message and exit
  -i INDIRNAME, --indir INDIRNAME
                        Input directory name were taxon_names, genes_taxonomy
                        and hmm_hits .txt files are found (required)
  -t TAXALEVEL, --taxa-level TAXALEVEL
                        Pick taxa level you want reported
  --print-taxa          Print taxa with nifHDK genes at level chosen to text
                        file.
  --get-taxa-contigs    Get contigs with nifHDK genes to text file.
  -g GENE, --gene GENE  Name of gene to search in dataframe for. You can pass
                        multiple genes by separating them by a pipe. Ex:
                        'geneA|geneB'
  --get-gene-contigs    Get contigs containing gene(s) of choice to text file.
 ```

To get the taxa for a particular gene call use the `anvi-get-taxonomy-for-genecall.py`


```
usage: anvi-get-taxonomy-for-genecall.py [-h] [-t TAXALEVEL] [--print-gene] -i INDIRNAME [-o OUTFILE] [-g GENE]

Parses anvi-export-table --table files

optional arguments:
  -h, --help            show this help message and exit
  -t TAXALEVEL, --taxa-level TAXALEVEL
                        Pick taxa level you want reported
  --print-gene          Print contigs with your choose combination of source and gene call(s).
  -i INDIRNAME, --indir INDIRNAME
                        Input directory name were taxon_names, genes_taxonomy and hmm_hits .txt files are found (required)
  -o OUTFILE, -outfile OUTFILE
                        Name of file to write output too.
  -g GENE, --gene-call GENE
                        Name of gene to search in dataframe for. You can pass multiple genes by separating them by a pipe. Ex: 'geneA|geneB'
```
