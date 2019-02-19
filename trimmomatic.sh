#!/bin/sh -l

#PBS -N trimm_1
#PBS -q christ99
#PBS -l nodes=1:ppn=4
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M 

module load trimmomatic
module load java

cd $PBS_O_WORKDIR


trimmomatic PE \
016401_Pool-1_ATCACG_R1.fastq \
016401_Pool-1_ATCACG_R2.fastq \
016401_Pool-1_ATCACG_R1_filtered.fastq \
016401_Pool-1_ATCACG_R1_unpaired.filtered.fastq \
016401_Pool-1_ATCACG_R2_filtered \
016401_Pool-1_ATCACG_R1_unpaired.filtered.fastq \
LEADING:21 TRAILING:20 MINLEN:30 ILLUMINACLIP:trimmomatic_illumina.fa:2:40:10


