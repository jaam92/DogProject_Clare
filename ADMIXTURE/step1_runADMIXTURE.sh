#!/bin/bash
#$ -cwd
#$ -V
#$ -N ADMIXTURE
#$ -l highp,h_data=4G,time=05:00:00
##$ -pe shared 4
#$ -M eplau
#$ -m bea

BED_FILE="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/PLINKFiles/prunedData_masterFile_allChroms.bed"
OUTDIR="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ADMIXTURE/outputFiles"
mkdir -p $OUTDIR

for K in {2..10} 
do

cd $OUTDIR 
/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/admixture_linux-1.3.0/admixture --cv $BED_FILE $K -j4 | tee log${K}.out
 
done

sleep 500
