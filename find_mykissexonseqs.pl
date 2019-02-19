#!/usr/bin/perl -w
#	A script to create a samtools script that finds sequence in reference genome.
#
#	Usage:
#	mykissexonsseqs.pl input.txt > outpuforjobsubmission.txt
#
#	October 2017 	Janna Willoughby
#	
#-----------------------------------------------------------------------------------------

use strict;

#set up variables
my $buffer = 51;
my $rm_out = $ARGV[0];
my $linenum = 0;
#set up file to read in
open(IN,"<$rm_out") || die ("Error opening $rm_out $!");

#print top lines for file 
print "#CDS records for $rm_out\n"; 

my $line = <IN>; 
while ($line = <>) {
	#count lines
	$linenum++;
	
	if($linenum<=1){
		next;
	}else{
		#prep data, extract chromosome type
		chomp $line;	
		my @field = split "\t", $line;
	
		#print output
		if($field[1]>$field[2]){
			my $start = $field[2] - $buffer;
			my $end   = $field[1] + $buffer;
			print "samtools faidx ../../lien_mykiss/omyV6/omyV6Chr.fasta $field[0]:$start-$end >> mykiss_exonsBUFFER$buffer.fa\n";
		}
		if($field[2]>$field[1]){
			my $start = $field[1] - $buffer;
			my $end   = $field[2] + $buffer;	
			print "samtools faidx ../../lien_mykiss/omyV6/omyV6Chr.fasta $field[0]:$start-$end >> mykiss_exonsBUFFER$buffer.fa\n";
		}
	}
}

;