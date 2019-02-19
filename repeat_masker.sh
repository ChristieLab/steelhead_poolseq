#!/bin/bash
#PBS -N repeat_masker
#PBS -q christ99
#PBS -l walltime=336:00:00
#PBS -l nodes=1:ppn=10
#PBS -l naccesspolicy=shared
#PBS -m ae
#PBS -M xx@purdue.edu

cd $PBS_O_WORKDIR

module load bioinfo
module load RepeatMasker/4.0.7

RepeatMasker omyV6Chr.fasta -species salmonidae -s -no_is -cutoff 255 -frag 20000 -pa 10