#!/bin/bash
#$ -cwd
#$ -V
#$ -N pwFst
#$ -l h_data=500M,time=05:00:00
#$ -M eplau
#$ -m bea

#load vcftools
. /u/local/Modules/default/init/modules.sh
module load vcftools

#Comput Pairwise Fst
vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/IR_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/AW_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/IR_vs_AW

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/IR_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/EW_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/IR_vs_EW

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/IR_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/BC_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/IR_vs_BC

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/IR_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/TM_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/IR_vs_TM

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/IR_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/PG_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/IR_vs_PG

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/IR_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/LB_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/IR_vs_LB

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/AW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/EW_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/AW_vs_EW

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/AW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/BC_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/AW_vs_BC

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/AW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/TM_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/AW_vs_TM

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/AW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/PG_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/AW_vs_PG

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/AW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/LB_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/AW_vs_LB

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/EW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/BC_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/EW_vs_BC

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/EW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/TM_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/EW_vs_TM

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/EW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/PG_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/EW_vs_PG

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/EW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/LB_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/EW_vs_LB

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/BC_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/TM_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/BC_vs_TM

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/BC_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/PG_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/BC_vs_PG

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/BC_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/LB_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/BC_vs_LB

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/TM_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/PG_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/TM_vs_PG

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/TM_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/LB_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/TM_vs_LB

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/PG_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/LB_indivList.txt  --out /u/flashscratch/j/jmooney3/pairwiseFst/PG_vs_LB
