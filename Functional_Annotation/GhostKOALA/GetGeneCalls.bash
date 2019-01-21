#!/bin/bash
##
#SBATCH --mem=20G
#SBATCH --time=3-18:0:0
#SBATCH -n 3
#SBATCH -p production
#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/GhostKOALA/Error_Out_Files/Get_GhostK.out
#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/GhostKOALA/Error_Out_Files/Get_GhostK.err

aklog
source activate anvio5.2

#get gene calls to run GhostKOALA
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.2/
time anvi-get-sequences-for-gene-calls --get-aa-sequences -c fixed.contigsV5_1000D.db -o ../GhostKOALA/protein-sequences.fa
echo "Done getting AA calls"
#adding genecall_ to the end of reach gene call so that the GhoastKOALA server isn't mad at us
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/GhostKOALA/
sed -i 's/>/>genecall_/g' protein-sequences.fa
echo "Done Changing gene call names"
