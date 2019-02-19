#!/usr/bin/perl -w
#
#	A script to combine alleles at a single location to one line, using polimaps output.
#
#	January 2018, Janna Willoughby
#	
#	note:see associated bash script that filters files by chromosome
#		
#-----------------------------------------------------------------------------------------

use strict;
use List::Util qw[min max];

#set up data and variables
my $all = $ARGV[0];

#initialize variables
my @all_locs;
my $prepos = 0;
my $chromosome;
my $allfreq = 0.05;

my $Afreq1   = 0;
my $Afreq2   = 0;
my $Afreq3   = 0;
my $ATdepth1 = 0;
my $ATdepth2 = 0;
my $ATdepth3 = 0;
my $Adepth1  = 0;
my $Adepth2  = 0;
my $Adepth3  = 0;

my $Tfreq1   = 0;
my $Tfreq2   = 0;
my $Tfreq3   = 0;
my $TTdepth1 = 0;
my $TTdepth2 = 0;
my $TTdepth3 = 0;
my $Tdepth1  = 0;
my $Tdepth2  = 0;
my $Tdepth3  = 0;

my $Cfreq1   = 0;
my $Cfreq2   = 0;
my $Cfreq3   = 0;
my $CTdepth1 = 0;
my $CTdepth2 = 0;
my $CTdepth3 = 0;
my $Cdepth1  = 0;
my $Cdepth2  = 0;
my $Cdepth3  = 0;

my $Gfreq1   = 0;
my $Gfreq2   = 0;
my $Gfreq3   = 0;
my $GTdepth1 = 0;
my $GTdepth2 = 0;
my $GTdepth3 = 0;
my $Gdepth1  = 0;
my $Gdepth2  = 0;
my $Gdepth3  = 0;

#print header for output file
#print "chromosome\tlocation\tA1\tT1\tC1\tG1\tA2\tT2\tC2\tG2\tA3\tT3\tC3\tG3\tdepthA1\tdepthT1\tdepthC1\tdepthG1\tdepthA2\tdepthT2\tdepthC2\tdepthG2\tdepthA3\tdepthT3\tdepthC3\tdepthG3\trefallele1\trefallele2\trefallele3\trefdepth1\trefdepth2\trefdepth3\ttotdepth1\ttotdepth2\ttotdepth3\n";

####reformat allele frequency information based on Fst locations####
open(ALL,"<$all") || die ("Error opening $all $!");   
my $line = <ALL>; 
while ($line = <>) {
	chomp $line;	
	my @cur_line       = split "\t", $line;
	my $cur_alleleinfo = $cur_line[0];
	my @cur_sepinfo    = split "_", $cur_alleleinfo;
	my $cur_location   = $cur_sepinfo[1];
	my $cur_nucl       = $cur_sepinfo[2];
	
	#check depth in each population, skip if < 20 reads
	if($cur_line[2] < 20 || $cur_line[4] < 20 || $cur_line[6] < 20){
		$prepos = $cur_location;
		next;
	}
	
	#set values in first line
	if($prepos==0){
		$prepos     = $cur_location;
		$chromosome = $cur_sepinfo[0];
	}
	
	#if allele location is larger than previous, process saved data and reset variables
	if($prepos < $cur_location){
		#check sums
		my $pop1 = $Afreq1 + $Tfreq1 + $Cfreq1 + $Gfreq1;
		my $pop2 = $Afreq2 + $Tfreq2 + $Cfreq2 + $Gfreq2;
		my $pop3 = $Afreq3 + $Tfreq3 + $Cfreq3 + $Gfreq3;
		
		if($pop1 > 0 || $pop2 > 0 || $pop3 > 0){
			###process saved data###
			#record total depth
			my $totdepth1  = max($ATdepth1 , $TTdepth1 , $CTdepth1 , $GTdepth1 );
			my $totdepth2  = max($ATdepth2 , $TTdepth2 , $CTdepth2 , $GTdepth2 );
			my $totdepth3  = max($ATdepth3 , $TTdepth3 , $CTdepth3 , $GTdepth3 );
			
			#calculate reference depth
			my $refdepth1  = $totdepth1 - ($Adepth1 + $Tdepth1 + $Cdepth1 + $Gdepth1);
			my $refdepth2  = $totdepth2 - ($Adepth2 + $Tdepth2 + $Cdepth2 + $Gdepth2);
			my $refdepth3  = $totdepth3 - ($Adepth3 + $Tdepth3 + $Cdepth3 + $Gdepth3);
			 
			#calculate reference frequency
			my $refallele1 = 1 - ($Afreq1 + $Tfreq1 + $Cfreq1 + $Gfreq1);
			my $refallele2 = 1 - ($Afreq2 + $Tfreq2 + $Cfreq2 + $Gfreq2);
			my $refallele3 = 1 - ($Afreq3 + $Tfreq3 + $Cfreq3 + $Gfreq3);
			
			#print output to file
			print "$chromosome\t$prepos\t$Afreq1\t$Tfreq1\t$Cfreq1\t$Gfreq1\t$Afreq2\t$Tfreq2\t$Cfreq2\t$Gfreq2\t$Afreq3\t$Tfreq3\t$Cfreq3\t$Gfreq3\t$Adepth1\t$Tdepth1\t$Cdepth1\t$Gdepth1\t$Adepth2\t$Tdepth2\t$Cdepth2\t$Gdepth2\t$Adepth3\t$Tdepth3\t$Cdepth3\t$Gdepth3\t$refallele1\t$refallele2\t$refallele3\t$refdepth1\t$refdepth2\t$refdepth3\t$totdepth1\t$totdepth2\t$totdepth3\n";
		}
		
		###reset variables###
		$Afreq1 = $Afreq2 = $Afreq3 = $Adepth1 = $Adepth2 = $Adepth3 = 0;
		$Tfreq1 = $Tfreq2 = $Tfreq3 = $Tdepth1 = $Tdepth2 = $Tdepth3 = 0;
		$Cfreq1 = $Cfreq2 = $Cfreq3 = $Cdepth1 = $Cdepth2 = $Cdepth3 = 0;
		$Gfreq1 = $Gfreq2 = $Gfreq3 = $Gdepth1 = $Gdepth2 = $Gdepth3 = 0;
	}
	
	#save information for future processing
	$prepos = $cur_location;
	$chromosome = $cur_sepinfo[0];
	
	if("A" =~ /$cur_nucl/){
		if($cur_line[1] > $allfreq){
 			$Afreq1   = $cur_line[1];
 			$Adepth1  = $cur_line[2] * $cur_line[1];
 		}
 		if($cur_line[3] > $allfreq){
 			$Afreq2   = $cur_line[3];
 			$Adepth2  = $cur_line[4] * $cur_line[3];
 		}
 		if($cur_line[5] > $allfreq){
 			$Afreq3   = $cur_line[5];
 			$Adepth3  = $cur_line[6] * $cur_line[5];
 		}
		$ATdepth1 = $cur_line[2];
		$ATdepth2 = $cur_line[4];
		$ATdepth3 = $cur_line[6];
		next;
	}
	if("T" =~ /$cur_nucl/){
		if($cur_line[1] > $allfreq){
 			$Tfreq1   = $cur_line[1];
 			$Tdepth1  = $cur_line[2] * $cur_line[1];
 		}
 		if($cur_line[3] > $allfreq){
 			$Tfreq2   = $cur_line[3];
 			$Tdepth2  = $cur_line[4] * $cur_line[3];
 		}
 		if($cur_line[5] > $allfreq){
 			$Tfreq3   = $cur_line[5];
 			$Tdepth3  = $cur_line[6] * $cur_line[5];
 		}
		$TTdepth1 = $cur_line[2];
		$TTdepth2 = $cur_line[4];
		$TTdepth3 = $cur_line[6];
		next;
	}
	if("C" =~ /$cur_nucl/){
		if($cur_line[1] > $allfreq){
 			$Cfreq1   = $cur_line[1];
 			$Cdepth1  = $cur_line[2] * $cur_line[1];
 		}
 		if($cur_line[3] > $allfreq){
 			$Cfreq2   = $cur_line[3];
 			$Cdepth2  = $cur_line[4] * $cur_line[3];
 		}
 		if($cur_line[5] > $allfreq){
 			$Cfreq3   = $cur_line[5];
 			$Cdepth3  = $cur_line[6] * $cur_line[5];
 		}
		$CTdepth1 = $cur_line[2];
		$CTdepth2 = $cur_line[4];
		$CTdepth3 = $cur_line[6];
		next;
	}
	if("G" =~ /$cur_nucl/){
		if($cur_line[1] > $allfreq){
 			$Gfreq1   = $cur_line[1];
 			$Gdepth1  = $cur_line[2] * $cur_line[1];
 		}
 		if($cur_line[3] > $allfreq){
 			$Gfreq2   = $cur_line[3];
 			$Gdepth2  = $cur_line[4] * $cur_line[3];
 		}
 		if($cur_line[5] > $allfreq){
 			$Gfreq3   = $cur_line[5];
 			$Gdepth3  = $cur_line[6] * $cur_line[5];
 		}
		$GTdepth1 = $cur_line[2];
		$GTdepth2 = $cur_line[4];
		$GTdepth3 = $cur_line[6];
		next;
	}
}
close ALL;
