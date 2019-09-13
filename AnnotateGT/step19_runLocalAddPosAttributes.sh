
#Takes about an hour to run
#If send to cluster you get a memory error

#Load python
. /u/local/Modules/default/init/modules.sh
module load python 

#set variables
for CHROM in {1..17}
do

for i in {AW10,AW11,AW12,AW13,AW15,AW16,AW17,AW18,AW19,AW20,AW21,AW22,AW23,AW24,BC1,BC3,BC4,BC5,BC6,BC7,BC9,EW1,EW10,EW2,EW3,EW4,EW5,EW6,EW7,EW8,EW9,IR1,IR10,IR2,IR3,IR4,IR5,IR6,IR7,IR8,IR9,LB1,LB11,LB2,LB3,LB4,LB5,LB6,LB7,LB8,LB9,PG1,PG10,PG11,PG12,PG13,PG14,PG16,PG2,PG3,PG4,PG5,PG6,PG7,PG8,PG9,TM1,TM10,TM2,TM3,TM5,TM6,TM7,TM8,TM9}
do

#give files
inputAnnotFile='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT/annotatedGTwithVEP_'$i'_Chr'$CHROM'.txt'
outputFile1='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT/annotatedGTwithVEP_'$i'_Chr'$CHROM'_addROHAnnot.txt'
outputFile2='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT/annotatedGTwithVEP_'$i'_Chr'$CHROM'_addROHAnnot_GERPScore.txt'

ROHrange1='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/vcfToolsROH/AnnotateROH/AnnotatedRange/SplitByIndivAndChrom/FinalROH_'$i'_chr'$CHROM'_Ranges.txt'
ROHrange2='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/vcfToolsROH/AnnotateROH/AnnotatedRange/SplitByIndivAndChrom_1MbGreater/FinalROH_'$i'_chr'$CHROM'_1MbRanges.txt'
GERPscores='/u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/CanFam3.1_gerpScoresDogProjClare/reformatForAnnotation_chr'$CHROM'.txt'

#Step 1: Intersect ROH Ranges
#python addPositionAttributes.py ${inputAnnotFile} 2 ${outputFile1} ranges ${ROHrange1} withinROH_10Kb ${ROHrange2} withinROH_1Mb  

#Step 2: Add GERP Score
python addPositionAttributes.py ${outputFile1} 2 ${outputFile2} positionAnnotation ${GERPscores} GERPScore

done
done


