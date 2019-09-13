#!/bin/bash
#$ -cwd 
#$ -V 
#$ -N copyVCFs
#$ -l time=08:00:00,h_data=500M 
#$ -M eplau 
#$ -m bea


for i in {AW,EW,IR,GS,PG,TM} 
do 

cp /u/home/c/cmarsden/nobackup-kirk/clare/analyses/FINAL_VCFS/final_filtered_vcf_allsites_6files/Ten"$i"/6_bespoke* /u/flashscratch/j/jmooney3/DogsROH/vcf_Ten"$i"/

done
