#!/bin/bash
#$ -l h_data=6G,h_rt=23:00:00
#$ -cwd
#$ -A jmooney3
#$ -N EWsims_3epoch
#$ -e /u/flashscratch/j/jmooney3
#$ -o /u/flashscratch/j/jmooney3
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load R/3.4.2
module load python

for i in {1..1000}
do

RandomNumber=$(( ( $SGE_TASK_ID - 1 ) * 1000 + $i ))

### Ethiopian Wolf Model modified from San Nicolas Model 5
python jazMod_RandomNumbersEW4_3epoch.py $SGE_TASK_ID ../InputData/EWSamps_allChroms_NeutralRegions_het_10000win_10000step_ABCinputFile.txt $RandomNumber 2000

###Grab new output file for each iteration through the randomnumbers script
File=$SCRATCH/"SimulationsAndParamsABC_3epoch/RandomNumbersEWsamps_"$SGE_TASK_ID".txt"
Output="HetEWsamps_"$SGE_TASK_ID".txt"

NumberOfRegions=$(  wc -l $File | awk '{print $1}' )

####Run the simulations using ms
/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msdir/ms 8 $NumberOfRegions -t tbs -r tbs tbs -eN tbs tbs -eN tbs tbs < $File | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msdir/sample_stats | awk '{print $2"\t"$4}' > msOutput_3epoch/$Output

####Compute the distance statistic for each iteration of simulations
Rscript --vanilla ComputeJointScore_ABC_3epoch_10KbWin.R --LengthThreshold 2000 --SGETaskID $SGE_TASK_ID --EmpInfile ../InputData/EWSamps_allChroms_NeutralRegions_het_10000win_10000step_ABCinputFile.txt  --SimInfile msOutput_3epoch/$Output --outFilePath /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/scripts_3epoch/ResultsABC_3epoch/ --paramFile $SCRATCH/SimulationsAndParamsABC_3epoch/OutputParams"$SGE_TASK_ID".txt


done

sleep 600
