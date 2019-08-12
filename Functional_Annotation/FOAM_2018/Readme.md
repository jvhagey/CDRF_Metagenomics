These are the instructions for how to take the full FOAM database of hmms and make a folder "Foam_Nitro" that you can point Anvi'o to in the anvi-run-hmm script  

1. FOAM-onto_rel1.tsv is the original file that comes with the foam database.  
2. Nitrogen_cycle_KO.txt is the list of KOs that are part of the nitrogen cycle from FOAM-onto_rel1.tsv  
3. To run anvi-run-hmms Anvi'o requires gene.txt, kind.txt, reference.txt, target.txt, noise_cutoff_terms.txt and genes.hmm.gz files. The FOAM_2018 folder has all the hmms from the FOAM database while the Foam_nitro folder has only a subset of the foam hmms that have ontology to the nitrogen cycle.  
4. To make the info in Foam_nitro folder run `grep -A1 "NAME" FOAM-hmm_rel1a.hmm > Acc_num.txt` in the terminal. You get the .hmm file from unzipping FOAM-hmm_rel1a.hmm.gz file.  
5. FOAM-hmm_rel1a.hmm has duplicate names :( like:   
```
NAME  KO:K00372_1.7.99.4
ACC   HMMsoil137302
NAME  KO:K00372_1.7.99.4
ACC   HMMsoil137304
NAME  KO:K00372_1.7.99.4
ACC   HMMsoil137400
```
Unfortunately, Anvi'o won't like this so you will need to run ```Edit_hmm.py -i ./input directory/ -o outputfile -i inputfile```  
6. Next, run the FOAM_to_anvio.py script to do a little parsing. The script should be run as `python FOAM_to_anvio.py -i ./input directory/`, the script will look for a file called Acc_num.txt in this folder. This generates the gene.txt file you will need.    
7. Nitro_KOs.txt is a list of unique KOs from Nitrogen_cycle_KO.tsv that we want to get out the FOAM database.  
8. Now run the hmmfetch.bash script that will allow you to get only the hmms you want from the FOAM-hmm_rel1a.hmm file.   
