#!/usr/bin/perl -w
#	A script to generate TE count data for windows of a set size with overlapping bins
#	of set size
#
#	Input = RepeatMasker .out file
#
#	July 2017			Avril Harder
#
#	modified October 2017 Janna Willoughby
#	Now parses more categories of output
#	note:remove top 3 lines (headers) before running!
#		example - tail -n +4 file.out > noheaderfile.out
#-----------------------------------------------------------------------------------------

use strict;

#set window size
my $window = 100000;

#set up variables
my $rm_out=$ARGV[0];
my @field;
my @class;

#these variables control window/position counting
my $low_lim = 1;
my $up_lim = $low_lim + $window;
my $win_num = 1;
my $s_pos;
my $chrom;

#counter for summarizing repeat types
my $sine_ct=0;	#SINE
my $line_ct=0;	#LINE (not Penelope)
my $pen_ct=0;	#LINE/Penelope
my $ltr_ct=0;	#LTR
my $dna_ct=0;	#DNA
my $srna_ct=0;	#snRNA
my $sat_ct=0;	#Satellite
my $simp_ct=0;	#Simple_repeat
my $low_ct=0;	#Low_complexity

#set up file to read in
open(IN,"<$rm_out") || die ("Error opening $rm_out $!");

#print header for output file
print "chromosome\tlocation\tSINEs\tLINES\tPenelopes\tLTRs\tDNAs\tsmallRNA\tsatellites\tsimple_reps\tlow_complex\n";

my $line = <IN>; 
while ($line = <>) {
	#print "$line\n";
 	chomp $line;	
	my @field = split " ", $line;
	my @class = split "\/", $field[10];
	
	# when a new chromosome is reached, start counters again
	if ($field[5] < $low_lim) {					
		$s_pos = 1;                
		$low_lim = 1;
		$up_lim = $low_lim + $window;
	
		#reset counts 
		$sine_ct=0;
		$line_ct=0;
		$pen_ct=0;
		$ltr_ct=0;
		$dna_ct=0;
		$srna_ct=0;
		$sat_ct=0;
		$simp_ct=0;
		$low_ct=0;
	}

	#look for window edge
	if ($field[5] >= $up_lim) {	
		#print output and increment variables
		$s_pos = $window * $win_num - $window + 1;
		print "$chrom\t$s_pos\t$sine_ct\t$line_ct\t$pen_ct\t$ltr_ct\t$dna_ct\t$srna_ct\t$sat_ct\t$simp_ct\t$low_ct\n";		
		$up_lim = $up_lim + $window;
		$low_lim = $low_lim + $window;
		$win_num++;
	
		#reset counts
		$sine_ct=0;
		$line_ct=0;
		$pen_ct=0;
		$ltr_ct=0;
		$dna_ct=0;
		$srna_ct=0;
		$sat_ct=0;
		$simp_ct=0;
		$low_ct=0;
	}

	#add occurence of each element type
	if ($field[5] >= $low_lim && $field[5] < $up_lim) {  
		if ($class[0] eq "SINE"){
			$sine_ct++;
		}
		if ($class[0] eq "LINE"){
			if ($class[1] eq "Penelope"){
				$pen_ct++;
			}else{
				$line_ct++;
			}
		}	
		if ($class[0] eq "LTR"){
			$ltr_ct++;
		}	
		if ($class[0] eq "DNA"){
			$dna_ct++;
		}
		if ($class[0] eq "snRNA"){
			$srna_ct++;
		}
		if ($class[0] eq "Satellite"){
			$sat_ct++;
		}
		if ($class[0] eq "Simple_repeat"){
			$simp_ct++;
		}
		if ($class[0] eq "Low_complexity"){
			$low_ct++;
		}
		$chrom = $field[4];	
	}	
}

;