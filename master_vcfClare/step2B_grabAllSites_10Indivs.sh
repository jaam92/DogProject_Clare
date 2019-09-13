#!/bin/bash
#$ -cwd
#$ -V
#$ -N pullSites10
#$ -l h_data=500M,time=12:00:00
#$ -M eplau
#$ -m bea

for i in {TenEW,TenIR,TenLB,TenTM,TenBC}
do 

for f in {1..38}
do

zcat /u/home/c/cmarsden/project-kirk-bigdata/FINAL_VCFS_May2018_INCrels/final_filtered_vcf_allsites_6files/"$i"/6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_"$i"_jointcalled_chr"$f".vcf.gz| grep -v "#" | awk '{print $2}' > /u/flashscratch/j/jmooney3/allSitesPerSpecies/"$i"_chr"$f"_allSitesPosOnly.txt

done
done

