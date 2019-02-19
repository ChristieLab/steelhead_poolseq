#!/bin/bash

#PBS -N PCoA
#PBS -q beagle
#PBS -l nodes=1:ppn=5,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m e
#PBS -M XX@purdue.edu


cd $PBS_O_WORKDIR

module load gcc/5.2.0 r/3.3.1
Rscript pcoa.R
