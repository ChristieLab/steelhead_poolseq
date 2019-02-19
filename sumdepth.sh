#!/bin/bash 

#PBS -N sum_depth3
#PBS -q beagle
#PBS -l nodes=1:ppn=2,naccesspolicy=shared
#PBS -l walltime=300:00:00
#PBS -m e
#PBS -M XX@purdue.edu

 
cd $PBS_O_WORKDIR

perl sumdepth_windows.pl depth_pop1AP.txt > sumdepth_1_AP.txt
perl sumdepth_windows.pl depth_pop2AP.txt > sumdepth_2_AP.txt
perl sumdepth_windows.pl depth_pop3AP.txt > sumdepth_3_AP.txt
