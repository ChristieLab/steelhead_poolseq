#!/bin/bash
 
#PBS -N picardsamtools_SS
#PBS -q christ99
#PBS -l nodes=1:ppn=1,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M XX@purdue.edu

cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load samtools/0.1.18

samtools mpileup -Df ../lien_mykiss/omyV6/omyV6Chr.fasta g1.sorted.all.newname.bam g2.sorted.all.newname.bam g3.sorted.all.newname.bam > Umykiss_g1g2g3.pileup

