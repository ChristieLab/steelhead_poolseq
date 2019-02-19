#!/bin/bash

#PBS -N picardsamtools_SS
#PBS -q christ99
#PBS -l nodes=1:ppn=6,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M XX@purdue.edu

cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load samtools/1.3.1

#convert sam to bam
samtools view -@ 10 -h -b -o g1.pool1.bam p1.sam
samtools view -@ 10 -h -b -o g1.pool2.bam p2.sam
samtools view -@ 10 -h -b -o g1.pool3.bam p3.sam
samtools view -@ 10 -h -b -o g1.pool4.bam p4.sam
samtools view -@ 10 -h -b -o g2.pool5.bam p5.sam
samtools view -@ 10 -h -b -o g2.pool6.bam p6.sam
samtools view -@ 10 -h -b -o g2.pool7.bam p7.sam
samtools view -@ 10 -h -b -o g2.pool8.bam p8.sam
samtools view -@ 10 -h -b -o g3.pool9.bam p9.sam
samtools view -@ 10 -h -b -o g3.pool10.bam p10.sam

#subsample reads
samtools view -@ 10 -s .99 -o g1.pool1.bam p1.bam
samtools view -@ 10 -o g1.pool2.bam p2.bam
samtools view -@ 10 -s .99 -o g1.pool3.bam p3.bam
samtools view -@ 10 -s .99 -o g1.pool4.bam p4.bam
samtools view -@ 10 -s .49 -o g2.pool5.bam p5.bam
samtools view -@ 10 -s .53 -o g2.pool6.bam p6.bam
samtools view -@ 10 -s .85 -o g2.pool7.bam p7.bam
samtools view -@ 10 -o g2.pool8.bam p8.bam
samtools view -@ 10 -o g3.pool9.bam p9.bam
samtools view -@ 10 -s .50 -o g3.pool10.bam p10.bam

#merge bams into groups
samtools merge -@ 10 -cp g1.all.bam g1.pool1.bam g1.pool2.bam g1.pool3.bam g1.pool4.bam
samtools merge -@ 10 -cp g2.all.bam g2.pool5.bam g2.pool6.bam g2.pool7.bam g2.pool8.bam
samtools merge -@ 10 -cp g3.all.bam g3.pool9.bam g3.pool10.bam 

#sort bams
samtools sort -@ 10 -o g1.sorted.all.bam g1.all.bam
samtools sort -@ 10 -o g2.sorted.all.bam g2.all.bam
samtools sort -@ 10 -o g3.sorted.all.bam g3.all.bam

#rename samples in bams
module load picard-tools/2.3.0
PicardCommandLine AddOrReplaceReadGroups I=g1.sorted.all.bam O=g1.sorted.all.newname.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=1 VALIDATION_STRINGENCY=SILENT
PicardCommandLine AddOrReplaceReadGroups I=g2.sorted.all.bam O=g2.sorted.all.newname.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=1 VALIDATION_STRINGENCY=SILENT
PicardCommandLine AddOrReplaceReadGroups I=g3.sorted.all.bam O=g3.sorted.all.newname.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=1 VALIDATION_STRINGENCY=SILENT
