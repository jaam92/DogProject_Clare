#!/bin/bash
#$ -l h_data=8G,h_rt=18:00:00
#$ -cwd
#$ -A jmooney3
#$ -N EWsims
#$ -e /u/flashscratch/j/jmooney3
#$ -o /u/flashscratch/j/jmooney3
#$ -m ea

. /u/local/Modules/default/init/modules.sh
module load R/3.4.2
module load python

for i in {1..1000}
do

RandomNumber=$(( ( $SGE_TASK_ID - 1 ) * 1000 + $i ))

### Ethiopian Wolf Model modified from San Nicolas Model 5
python jazMod_RandomNumbersEW4_5.py $SGE_TASK_ID /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt $RandomNumber 200

###Grab new output file for each iteration through the randomnumbers script
File=$SCRATCH/"SimulationsAndParamsABC/RandomNumbersEWsamps_"$SGE_TASK_ID".txt"
Output="HetEWsamps_"$SGE_TASK_ID".txt"

NumberOfRegions=$(  wc -l $File | awk '{print $1}' )

####Run the simulations using ms
/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msdir/ms 8 $NumberOfRegions -t tbs -r tbs tbs -eN tbs tbs  < $File | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msdir/sample_stats | awk '{print $2"\t"$4}' > msOutput_2epoch/$Output

####Compute the distance statistic for each iteration of simulations
Rscript --vanilla ComputeJointScore_ABC.R --LengthThreshold 200 --SGETaskID $SGE_TASK_ID --EmpInfile /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt  --SimInfile msOutput_2epoch/$Output --outFilePath /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC//u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/scripts_2epoch/ResultsABC/ --paramFile $SCRATCH/SimulationsAndParamsABC/OutputParams"$SGE_TASK_ID".txt


done

sleep 600
