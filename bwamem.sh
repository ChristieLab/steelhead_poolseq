#!/bin/bash

#PBS -N bwamem_ssalar
#PBS -q christ99
#PBS -l nodes=1:ppn=10,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M XX@purdue.edu

cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load bwa/0.7.10

#create index
bwa index -p mykiss_lien -a bwtsw ../lien_mykiss/omyV6/omyV6Chr.fasta &> bwaindex.log

#map reads
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016401_Pool-1_ATCACG_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016401_Pool-1_ATCACG_R2_filteredNseq.fastq > mykiss_Pool-1_ATCACG.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016402_Pool-2_CGATGT_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016402_Pool-2_CGATGT_R2_filteredNseq.fastq > mykiss_Pool-2_CGATGT.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016403_Pool-3_TTAGGC_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016403_Pool-3_TTAGGC_R2_filteredNseq.fastq > mykiss_Pool-3_TTAGGC.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016404_Pool-4_TGACCA_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016404_Pool-4_TGACCA_R2_filteredNseq.fastq > mykiss_Pool-4_TGACCA.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016405_Pool-5_ACAGTG_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016405_Pool-5_ACAGTG_R2_filteredNseq.fastq > mykiss_Pool-5_ACAGTG.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016406_Pool-6_GCCAAT_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016406_Pool-6_GCCAAT_R2_filteredNseq.fastq > mykiss_Pool-6_GCCAAT.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016407_Pool-7_CAGATC_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016407_Pool-7_CAGATC_R2_filteredNseq.fastq > mykiss_Pool-7_CAGATC.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016408_Pool-8_ACTTGA_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016408_Pool-8_ACTTGA_R2_filteredNseq.fastq > mykiss_Pool-8_ACTTGA.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016409_Pool-9_GATCAG_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016409_Pool-9_GATCAG_R2_filteredNseq.fastq > mykiss_Pool-9_GATCAG.sam
bwa mem -C -t 10 mykiss_lien /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016410_Pool-10_TAGCTT_R1_filteredNseq.fastq /scratch/snyder/j/jwillou/salmon_poolseq/data_programs/cleanedreads/renamed/016410_Pool-10_TAGCTT_R2_filteredNseq.fastq > mykiss_Pool-10_TAGCTT.sam
