#!/usr/bin/perl -w

#PBS -N separate_01
#PBS -q christ99
#PBS -l nodes=1:ppn=1,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M XX@purdue.edu

#separates single mapped reads from multiple mapped reads in sam file

chdir '/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/filteredreads';
#chdir '/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/ssalar/sambam/allmapped/mappedreads/perlmap';

my $file       = "../bwa/mykiss_Pool-1_ATCACG.sam";
my $outfile_h  = "pool1.header.sam";
my $outfile_u  = "pool1.unique.sam";
my $outfile_nu = "pool1.nonunique.sam";
my $countfile  = "pool1.numbermapped.txt";
open (INFILE,     "<", $file);
open (OUTFILE_H,  ">", $outfile_h);
open (OUTFILE_U,  ">", $outfile_u);
open (OUTFILE_NU, ">", $outfile_nu);
open (COUNTFILE,  ">", $countfile);


#read in file, line by line
while(my $line = <INFILE>){
	
	#split on tab
	my @sline   = split("\t", $line);
	
	#check to see if it is header line
	if($sline[0] =~ m/^@/){
		print OUTFILE_H $line;
		next;
	}
	
	#extract column of interest
	my $tomatch = $sline[15];
	
	#look for occurrence of NC in column (indicates additional mappings)
	my @matched = $tomatch =~ /omy/g;
	my $nmatch  = @matched;
	
	if($nmatch > 0){
		#if mapped multiple times, print to 'non unique' file
		print OUTFILE_NU $line;
		#count number of mapped times, add to file
		$nmatch++;
		print COUNTFILE "$nmatch \n";
		
	}else{
		#if mapped once times, print to 'unique' file
		print OUTFILE_U $line;
		#add single mapping to the number of mapped locations file
		print COUNTFILE "1 \n";
	}
}

close INFILE;
close OUTFILE_H;
close OUTFILE_U;
close OUTFILE_NU;
close COUNTFILE;

exit 0;
