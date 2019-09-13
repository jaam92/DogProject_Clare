#!/bin/bash
#$ -cwd
#$ -V
#$ -N Perms
#$ -l highp,h_data=5G,time=10:00:00 
#$ -pe shared 4
#$ -M eplau
#$ -m bea

#Generate Final counts of annotations
. /u/local/Modules/default/init/modules.sh
module load R/3.4.2

#Run scripts (script is set to be split across 4 cores)
#Rscript newVersion_ChristianPerm_1Mb.R
#Rscript newVersion_ChristianPerm.R
#Rscript newVersion_ChristianPerm_btwnTenKbandOneMb.R 


#Rscript newVersion_ChristianPerm_1Mb_GERP.R


