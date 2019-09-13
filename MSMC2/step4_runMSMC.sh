#!/bin/bash
#$ -cwd 
#$ -V 
#$ -N msmc2
#$ -l highp,time=15:00:00,h_data=4G
##$ -pe 4 
#$ -M eplau 
#$ -m bea

for j in {EW2,EW3,EW9,EW10,AW11,AW20,IR5,IR10}
do
/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msmc2-master/./msmc2 -i 25 -t 4 -o /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/MSMC2/outputMSMC/msmc2_output_iter50_"$j" "$j"_chr1_msmcInput.txt "$j"_chr2_msmcInput.txt "$j"_chr3_msmcInput.txt "$j"_chr4_msmcInput.txt "$j"_chr5_msmcInput.txt "$j"_chr6_msmcInput.txt "$j"_chr7_msmcInput.txt "$j"_chr8_msmcInput.txt "$j"_chr9_msmcInput.txt "$j"_chr10_msmcInput.txt "$j"_chr11_msmcInput.txt "$j"_chr12_msmcInput.txt "$j"_chr13_msmcInput.txt "$j"_chr14_msmcInput.txt "$j"_chr15_msmcInput.txt "$j"_chr16_msmcInput.txt "$j"_chr17_msmcInput.txt "$j"_chr18_msmcInput.txt "$j"_chr19_msmcInput.txt "$j"_chr20_msmcInput.txt "$j"_chr21_msmcInput.txt "$j"_chr22_msmcInput.txt

/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msmc2-master/./msmc2 -i 50 -t 4 -o /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/MSMC2/outputMSMC/msmc2_output_iter100_"$j" "$j"_chr1_msmcInput.txt "$j"_chr2_msmcInput.txt "$j"_chr3_msmcInput.txt "$j"_chr4_msmcInput.txt "$j"_chr5_msmcInput.txt "$j"_chr6_msmcInput.txt "$j"_chr7_msmcInput.txt "$j"_chr8_msmcInput.txt "$j"_chr9_msmcInput.txt "$j"_chr10_msmcInput.txt "$j"_chr11_msmcInput.txt "$j"_chr12_msmcInput.txt "$j"_chr13_msmcInput.txt "$j"_chr14_msmcInput.txt "$j"_chr15_msmcInput.txt "$j"_chr16_msmcInput.txt "$j"_chr17_msmcInput.txt "$j"_chr18_msmcInput.txt "$j"_chr19_msmcInput.txt "$j"_chr20_msmcInput.txt "$j"_chr21_msmcInput.txt "$j"_chr22_msmcInput.txt

done


