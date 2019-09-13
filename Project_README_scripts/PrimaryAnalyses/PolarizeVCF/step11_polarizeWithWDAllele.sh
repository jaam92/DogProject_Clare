#!/bin/bash
#$ -cwd
#$ -V
#$ -N RecodeVCF
#$ -l time=00:40:00,h_data=4G
#$ -M eplau
#$ -m bea

###Each step takes about 20 minutes to run locally

script_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RecodedVCF'
refAllele_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/WildDog'
in_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RecodedVCF/GenotypeOnly'
out_dir_WD='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RecodedVCF/WDAllele' 
out_dir_GT='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RecodedVCF'

for CHROM in {1..38}
do

MasterFile_vcfin=${in_dir}/'reformatted_masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.vcf'

WildDog_allele=${refAllele_dir}/'WD_Alleles_inMasterFile_Chr'${CHROM}'.txt'

WildDogRecode_vcfout=${out_dir_WD}/'WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.vcf'

summaryFile=${out_dir_WD}/'summaryRecode_WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.out'

genoStateRecode_vcfout=${out_dir_GT}/'AssignGT_WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.vcf'

#Generate True ROH
. /u/local/Modules/default/init/modules.sh
module load python 

#Recode with Wild Dog Allele
#python ${script_dir}/Recode_withWildDog_ref_final.py ${WildDog_allele} ${MasterFile_vcfin} ${WildDogRecode_vcfout} > ${summaryFile}

#Write in the new genotype
python ${script_dir}/assign_genotype_state.py ${WildDogRecode_vcfout} ${genoStateRecode_vcfout}

done
