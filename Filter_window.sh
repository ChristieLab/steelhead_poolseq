#!/bin/bash 

#PBS -N filter_mykiss
#PBS -q christ99
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=300:00:00
#PBS -m e
#PBS -M jwillou@purdue.edu

 
cd $PBS_O_WORKDIR

module load gcc/5.2.0 r/3.3.1
Rscript Filter_window.R

mkdir plots
Rscript findandplot_outliers.R

