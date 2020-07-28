#!/bin/bash
##
#SBATCH --mem=5G
#SBATCH -p production
#SBATCH -n 3
#SBATCH --time=6-18:0:0
#SBATCH -o /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/Error_Out_Files/metabat_anvio_contigs_nif.out
#SBATCH -e /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/SourMash/Error_Out_Files/metabat_anvio_contigs_nif.err


module load seqtk/1.3

##first look for contigs have contain nifHDK in the same contig
#using seqk search fasta file for specific contigs.
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/SourMash/
time seqtk subseq ../../C5_V5.1/fixed.contigsV5_1000.fa Allmethods_HDK_contigs.txt > contigs_HDK_Allmethods.fa

grep -H -f contigs_HDK_Allmethods.fa ../../../MetaBat/C5_bin_1500/*.fa > bins_with_HDK_Allmethods.txt
#use this to get contig names alone
sed 's/:>c_[0-9]\+$//' bins_with_HDK_Allmethods.txt | sed 's/.*C5_bin_1500\///' | uniq > bins_with_HDK_Allmethods2.txt

mkdir MetaBat_bins_Allmethods

#copy the fasta files into new folder.
cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/MetaBat/C5_bin_1500/

cat ../../Anvio/C5_V5.5/SourMash/bins_with_HDK_Allmethods2.txt | xargs cp -t ../../Anvio/C5_V5.5/SourMash/MetaBat_bins_Allmethods/

##finding bins that have contigs that have each gene individually 
time seqtk subseq ../../C5_V5.1/fixed.contigsV5_1000.fa Allmethods_nifH_contigs.txt > contigs_nifH_Allmethods.fa
time seqtk subseq ../../C5_V5.1/fixed.contigsV5_1000.fa Allmethods_nifKD_contigs.txt > contigs_nifDK_Allmethods.fa
grep -H -f contigs_nifH_Allmethods.fa ../../../MetaBat/C5_bin_1500/*.fa > bins_with_nifH_Allmethods.txt
grep -H -f contigs_nifDK_Allmethods.fa ../../../MetaBat/C5_bin_1500/*.fa > bins_with_nifDK_Allmethods.txt
#use this to get contig names alone
sed 's/:>c_[0-9]\+$//' bins_with_nifH_Allmethods.txt | sed 's/.*C5_bin_1500\///' | uniq > bins_with_nifH_Allmethods2.txt
sed 's/:>c_[0-9]\+$//' bins_with_nifDK_Allmethods.txt | sed 's/.*C5_bin_1500\///' | uniq > bins_with_nifDK_Allmethods2.txt
#make directories
mkdir ../../Anvio/C5_V5.5/SourMash/MetaBat_bins_nifH_Allmethods
mkdir ../../Anvio/C5_V5.5/SourMash/MetaBat_bins_nifDK_Allmethods
#copy the fasta files into new folder.
cat ../../Anvio/C5_V5.5/SourMash/bins_with_nifH_Allmethods2.txt | xargs cp -t ../../Anvio/C5_V5.5/SourMash/MetaBat_bins_nifH_Allmethods/
cat ../../Anvio/C5_V5.5/SourMash/bins_with_nifDK_Allmethods2.txt | xargs cp -t ../../Anvio/C5_V5.5/SourMash/MetaBat_bins_nifDK_Allmethods/

cd /share/tearlab/Maga/Jill/CDRF_MetaGenome/Assembly_2018/Anvio/C5_V5.5/SourMash/
#removing some bins we don't want to look at
rm ./MetaBat_bins_nifDK_Allmethods/C5_bin_more.unbinned.fa
rm ./MetaBat_bins_nifH_Allmethods/C5_bin_more.unbinned.fa
#we alreadly looked at these bins so removing them too
rm ./MetaBat_bins_nifDK_Allmethods/C5_bin_more.205.fa
rm ./MetaBat_bins_nifH_Allmethods/C5_bin_more.205.fa
rm ./MetaBat_bins_nifDK_Allmethods/C5_bin_more.173.fa
rm ./MetaBat_bins_nifH_Allmethods/C5_bin_more.173.fa
rm ./MetaBat_bins_nifDK_Allmethods/C5_bin_more.57.fa
rm ./MetaBat_bins_nifH_Allmethods/C5_bin_more.57.fa
rm ./MetaBat_bins_nifDK_Allmethods/C5_bin_more.83.fa
rm ./MetaBat_bins_nifH_Allmethods/C5_bin_more.83.fa
#check for common bins
diff -srq MetaBat_bins_nifDK_Allmethods/ MetaBat_bins_nifH_Allmethods/ | grep identical > commmon_bins.txt
