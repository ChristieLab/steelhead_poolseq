#!/usr/bin/perl -w
#-----------------------------------------------------------------------------------------
#   calculate depth size from bedtools output
#   example bedtools:
#   bedtools genomecov -d -ibam test.sorted.bam -g omyV6Chr.fasta > depth_test.txt
#   usage: perl sumdepth_windows.pl depth_test.txt depth_windows.txt
#   October 2017, Janna Willoughby
#-----------------------------------------------------------------------------------------
#	Modified from 1start_rm_out_to_te_freqs.pl:
#   A script to generate TE count data for windows of a set size with overlapping bins
#	of set size
#	Input = RepeatMasker .out file
#	July 2017			Avril Harder
#-----------------------------------------------------------------------------------------

use strict;

#set window size
my $window = 100000;

#set up variables
my $rm_out=$ARGV[0];
my @field;

#these variables control window/position counting
my $low_lim = 1;                     #lower window position 
my $up_lim = $low_lim + $window;     #upper window position 
my $s_pos = 1;                       #current position

#counter for summing depths within windows
my $depth_ct = 0;

open (IN, "<$rm_out") || die ("Error opening $rm_out $!");
#open (OUT, '>depthwindows.txt' ) or die "cannot open output file";
#not smart enough to get the output file to work today, will need to specify out file command line

#print header for output file
print "chromosome\tlocation\tdepth\n";

#set infile, begin looping over lines
my $line = <IN>; 
while ($line = <>) {
 	chomp $line;	
	@field = split "\t", $line;
	
	# when a new chromosome is reached, start counters again
	if ($field[1] < $low_lim) {					
		$s_pos = 1;                
		$low_lim = 1;
		$up_lim = $low_lim + $window;
		$depth_ct = 0;
	}
	
	# sum depth with previous depth counts
	$depth_ct = $depth_ct + $field[2];
	
	# when move past current window boundaries ($field[1]), print totals begin move to the next window
	if ($field[1] >= $up_lim){
		#calculate average
		my $avgdepth = $depth_ct/$window;
		print "$field[0]\t$s_pos\t$avgdepth\n";		
		
		#set new window limits
		$up_lim  = $up_lim + $window;
		$low_lim = $low_lim + $window;
		
		#reset depth counter
		$depth_ct = 0;
		
		#set postion to next genome position
		$s_pos = $field[1]
	
	}
}
;