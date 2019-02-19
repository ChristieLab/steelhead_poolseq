#!/bin/bash

#PBS -N filterNs
#PBS -q christ99
#PBS -l nodes=1:ppn=1,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M X@purdue.edu

cd $PBS_O_WORKDIR

/usr/bin/perl -w FilterFastqPE.pl -f inputfilelist.txt -q 4 -g 75 -t 1 
