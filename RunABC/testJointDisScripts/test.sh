

#python jazMod_RandomNumbersEW4_5.py 38 ../InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt 37001 200
#wc -l RandomNumbersEWsamps_38.txt #use this number as second argument for ms
#/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msdir/ms 8 21724 -t tbs -r tbs tbs -eN tbs tbs < RandomNumbersEWsamps_38.txt | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msdir/sample_stats | awk '{print $2"\t"$4}' > HetEWsamps_38.txt

#run this step twice ( once for Diegos distance measure the second to normalize by empirical)
#python jazMod_CompareDistributions_UpperBinv3_PiandS.py ../InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt 200 HetEWsamps_38.txt 38 0 0 3 30 30 /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/jointDis/ /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/jointDis/ABCOutput_NoError_EWsamps_

#Run the new R script that is the equivalent of Diegos script
#. /u/local/Modules/default/init/modules.sh
#module load R/3.4.2
#Rscript --vanilla CompareBins.R --LengthThreshold 200 --SGETaskID 38 --EmpInfile ../InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt  --SimInfile HetEWsamps_38.txt --outFilePath /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/jointDis/ --paramFile OutputParams38.txt


#The Score for segsites should match across methods 
#I chose to compute score using diego's method so the score output from the R script matches the score from the first line of the output from the python script -> value is 22040 
