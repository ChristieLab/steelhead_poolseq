#!/bin/bash

#PBS -N allele_fstP2
#PBS -q christ99
#PBS -l nodes=1:ppn=7,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M XX@purdue.edu

cd $PBS_O_WORKDIR

perl AlleleFreqsFromPileupGeneralized.pl -g Umykiss_g1g2g3P2.pileup -d 20 -m 1
perl FstFromPooledFreqs.pl -f Allele_Frequencies_Umykiss_g1g2g3P2.pileup -d 20 -t 0.0001
