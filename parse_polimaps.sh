#!/bin/bash 

#PBS -N Parse_polimaps
#PBS -q darwin
#PBS -l nodes=1:ppn=5,naccesspolicy=shared
#PBS -l walltime=336:00:00
#PBS -m abe
#PBS -M XX@purdue.edu

cd $PBS_O_WORKDIR

#break files down by chromosome
cat chromosomes.txt | while read -r LINE
	do
	name=$(echo $LINE | awk '{print $1}')
	grep  $name ../../callSNPs/Allele_Frequencies_Umykiss_g1g2g3.pileup > Alleles_temp.txt
	perl parse_polimaps_combined.pl Alleles_temp.txt > $name	
done

#combine output into single file
cat omy01 omy02 omy03 omy04 omy05 omy06 omy07 omy08 omy09 omy10 omy11 omy12 omy13 omy14 omy15 omy16 omy17 omy18 omy19 omy20 omy21 omy22 omy23 omy24 omy25 omy26 omy27 omy28 omy29 > Allele_Frequencies_Umykiss_g1g2g3.pileup_parsed
