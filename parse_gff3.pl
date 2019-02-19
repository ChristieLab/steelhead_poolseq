#!/usr/bin/perl -w
#	A script to parse annotation (gff3) files and extract only CDS locations. Modified 
#	from Repeat Masker parsing script (Avril Harder, July 2017). 
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
my @field;
my $class;
my $linenum = 0;
#set up file to read in
open(IN,"<$rm_out") || die ("Error opening $rm_out $!");

#print header for output file 
print "#CDS records for $rm_out\n"; 

my $line = <IN>; 
while ($line = <>) {
 	#prep data, extract annotation type
 	chomp $line;	
	my @field = split "\t", $line;
	my $class = $field[2];
	
	#count lines
	$linenum++;
	
	#skip header line
	if($linenum eq 1){
		next;
	}else{
		if($class eq "CDS"){
			print "$field[0]\t$field[1]\t$field[2]\t$field[3]\t$field[4]\n";
		}
	}
}

;