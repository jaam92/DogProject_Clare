#$ -l h_rt=20:00:00
#$ -cwd
#$ -A jmooney3
#$ -N EWsims
#$ -e /u/flashscratch/j/jmooney3
#$ -o /u/flashscratch/j/jmooney3

#Generate the output file and write header line to it
#header line will change based on number of parameters being estimated and maximum number of heterozygote bins used
ABCOutputFile="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/ABCResults/ABCOutput_NoError_EWsamps_"$SGE_TASK_ID".txt"

touch $ABCOutputFile

for i in {1..1000}
do

RandomNumber=$(( ( $SGE_TASK_ID - 1 ) * 1000 + $i ))

### Ethiopian Wolf Model modified from San Nicolas Model 5
python jazMod_RandomNumbersEW4_5.py $SGE_TASK_ID InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt $RandomNumber 200

###Grab new output file for each iteration through the randomnumbers script
File=$SCRATCH/"SimulationsAndParamsABC/RandomNumbersEWsamps_"$SGE_TASK_ID".txt"
Output="HetEWsamps_"$SGE_TASK_ID".txt"

NumberOfRegions=$(  wc -l $File | awk '{print $1}' )

####Run the simulations using ms
/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msdir/ms 2 $NumberOfRegions -t tbs -r tbs tbs -eN tbs tbs  < $File | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msdir/sample_stats | awk '{print $4}' > msOutput/$Output

####Compute the distance statistic for each iteration of simulations
####maximum hets is set to 90 since we are merging 4 individuals data together 
python jazMod_CompareDistributions_UpperBinv3.py InputData/EW4_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt 200 msOutput/$Output $SGE_TASK_ID 0 0 3 25 25 /u/flashscratch/j/jmooney3/SimulationsAndParamsABC/ /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/ABCResults/ABCOutput_NoError_EWsamps_

done

sleep 300
