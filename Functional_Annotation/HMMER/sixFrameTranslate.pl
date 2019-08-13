#!/usr/bin/perl
#Thanks to Dr. Srijak Bhatnagar for this script!
use strict;
use Bio::SeqIO;
use Bio::SeqUtils;
unless(@ARGV){	die "\nUSAGE: $0 -i <input_file> -o <output_file>\n";	}
my %arg;
%arg = @ARGV;
my $seqIn = Bio::SeqIO->new(-file => $arg{'-i'});
my $seqOut = Bio::SeqIO->new(-file   => ">$arg{'-o'}", -format => 'fasta');
while(my $seq = $seqIn->next_seq() ) {
	my @seqs;
	@seqs = Bio::SeqUtils->translate_6frames($seq);
	print $seqOut->write_seq($seqs[0]);
	print $seqOut->write_seq($seqs[1]);
	print $seqOut->write_seq($seqs[2]);
	print $seqOut->write_seq($seqs[3]);
	print $seqOut->write_seq($seqs[4]);
	print $seqOut->write_seq($seqs[5]);
}
