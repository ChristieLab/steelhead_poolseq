#!/bin/bash

#PBS -N assemble
#PBS -q christ99
#PBS -l nodes=1:ppn=1,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M XX@purdue.edu

cd $PBS_O_WORKDIR
 
cat pool1.header.sam pool1.unique.sam > p1.sam
cat pool2.header.sam pool2.unique.sam > p2.sam
cat pool3.header.sam pool3.unique.sam > p3.sam
cat pool4.header.sam pool4.unique.sam > p4.sam
cat pool5.header.sam pool5.unique.sam > p5.sam
cat pool6.header.sam pool6.unique.sam > p6.sam
cat pool7.header.sam pool7.unique.sam > p7.sam
cat pool8.header.sam pool8.unique.sam > p8.sam
cat pool9.header.sam pool9.unique.sam > p9.sam
cat pool10.header.sam pool10.unique.sam > p10.sam
