#!/bin/bash
##
#SBATCH -p test
#SBATCH --mem=50G
#SBATCH --time=5-18:0:0
#SBATCH -n 5

module load hmmer/3.1b2

cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/FOAM_2018/
#Make list of names to extract, sort and return unique names
grep -f Nitro_KOs.txt FOAM-hmm_rel1a.hmm | sort -u > Nitro_KOs_Full.txt
#remove "NAME  " at the beginning of the line
time hmmfetch --index FOAM-hmm_rel1a.hmm
time hmmfetch -o Nitro_Foam.hmm -f FOAM-hmm_rel1a.hmm Nitro_KOs_Full.txt
