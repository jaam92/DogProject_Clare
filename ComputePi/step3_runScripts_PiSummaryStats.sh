#!/bin/bash
#$ -cwd
#$ -V
#$ -N summaryStats
#$ -l highp,h_data=6G,time=03:00:00 
#$ -pe shared 4
#$ -M eplau
#$ -m bea

#Generate Final counts of annotations
. /u/local/Modules/default/init/modules.sh
module load R/3.4.2

#Run scripts (script is set to be split across 4 cores)
Rscript SummarizePi_ROHvsnonROH_allPops.R

