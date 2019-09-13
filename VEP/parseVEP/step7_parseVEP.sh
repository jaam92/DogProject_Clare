#!/bin/bash
#$ -cwd
#$ -V
#$ -N parseVEP
#$ -l time=10:00:00,h_data=5G
#$ -M eplau
#$ -m bea

#Generate Final VEP annotations
. /u/local/Modules/default/init/modules.sh
module load R/3.4.2

Rscript parseVEPOutput_v2.R
