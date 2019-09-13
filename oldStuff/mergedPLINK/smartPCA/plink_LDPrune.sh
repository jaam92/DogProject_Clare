#!/bin/bash
#$ -cwd 
#$ -V 
#$ -N prune
#$ -l h_data=1G,time=00:20:00 
#$ -M eplau 
#$ -m bea

#can run locally do not need to qsub
#takes about 5mins to do each step

#Load Plink
source /u/local/Modules/default/init/bash
module load plink

work_dir='/u/flashscratch/j/jmooney3/DogsROH/mergedPLINK/smartPCA'

plink_in=${work_dir}/'allDogsMerged_smartPCA'
prune_in=${work_dir}/'keepSites_LDPrune'

plink_out=${work_dir}/'allDogsMerged_smartPCA_LDpruned'

#plink_Fsnp=${work_dir}/'allDogsMerged_smartPCA_LDpruned_Fsnp'

	#LD prune snps 
	#plink --bfile ${plink_in} --maf 0.05 --indep-pairwise 50 5 0.5 --dog --out ${prune_in} 


	#Remove sites make new file
	#plink --bfile ${plink_in} --extract keepSites_LDPrune.prune.in --make-bed --dog --out ${plink_out}

	#Calculate Fsnp
	#plink --bfile ${plink_out} --het --dog --out ${plink_Fsnp}

#sleep 300 

