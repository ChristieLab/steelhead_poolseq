#!/usr/bin/perl -w
#	A script to  create a samtools script that finds sequence in reference genome.
#
#	Usage:
#	parse_gff3.pl	input.gff3 > CDSinput.gff3.txt
#
#	October 2017 	Janna Willoughby
#	
#-----------------------------------------------------------------------------------------

use strict;

#set up variables
my $rm_out=$ARGV[0];
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
		my $chrom = $field[0];
		#print "$chrom\n";
	
		#replace chromosome designation in the longest and stupidest way possible (other than by hand)
		if($chrom =~ /ssa01/){$chrom='NC_027300_chr01';}
		if($chrom =~ /ssa02/){$chrom='NC_027301_chr02';}
		if($chrom =~ /ssa03/){$chrom='NC_027302_chr03';}
		if($chrom =~ /ssa04/){$chrom='NC_027303_chr04';}
		if($chrom =~ /ssa05/){$chrom='NC_027304_chr05';}
		if($chrom =~ /ssa06/){$chrom='NC_027305_chr06';}
		if($chrom =~ /ssa07/){$chrom='NC_027306_chr07';}
		if($chrom =~ /ssa08/){$chrom='NC_027307_chr08';}
		if($chrom =~ /ssa09/){$chrom='NC_027308_chr09';}
		if($chrom =~ /ssa10/){$chrom='NC_027309_chr10';}
		if($chrom =~ /ssa11/){$chrom='NC_027310_chr11';}
		if($chrom =~ /ssa12/){$chrom='NC_027311_chr12';}
		if($chrom =~ /ssa13/){$chrom='NC_027312_chr13';}
		if($chrom =~ /ssa14/){$chrom='NC_027313_chr14';}
		if($chrom =~ /ssa15/){$chrom='NC_027314_chr15';}
		if($chrom =~ /ssa16/){$chrom='NC_027315_chr16';}
		if($chrom =~ /ssa17/){$chrom='NC_027316_chr17';}
		if($chrom =~ /ssa18/){$chrom='NC_027317_chr18';}
		if($chrom =~ /ssa19/){$chrom='NC_027318_chr19';}
		if($chrom =~ /ssa20/){$chrom='NC_027319_chr20';}
		if($chrom =~ /ssa21/){$chrom='NC_027320_chr21';}
		if($chrom =~ /ssa22/){$chrom='NC_027321_chr22';}
		if($chrom =~ /ssa23/){$chrom='NC_027322_chr23';}
		if($chrom =~ /ssa24/){$chrom='NC_027323_chr24';}
		if($chrom =~ /ssa25/){$chrom='NC_027324_chr25';}
		if($chrom =~ /ssa26/){$chrom='NC_027325_chr26';}
		if($chrom =~ /ssa27/){$chrom='NC_027326_chr27';}
		if($chrom =~ /ssa28/){$chrom='NC_027327_chr28';}
		if($chrom =~ /ssa29/){$chrom='NC_027328_chr29';}
	
		#print output
		print "samtools faidx ../../ssalar/genome/ssalar_chr.fa $chrom:$field[3]-$field[4] >> ssalar_CDS.fa\n";
	}
}

;