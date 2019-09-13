import sys 
import random 
import math

#Input Arguments
for arg in sys.argv:
	print arg
Number = int(sys.argv[1])
LengthFile = sys.argv[2]
RandomSeed = int(sys.argv[3])
LengthThreshold = int(sys.argv[4])

#Initiate Counters and Lists and variables
random.seed( RandomSeed )
MutRate = 0.00000000496 #Tanya Autosomal Mutation Rate for Neutral cm >0.4
Reps = 0;
CalledSites = []

#Output Files
RandomNumbersFile = "/u/flashscratch/j/jmooney3/SimulationsAndParamsABC/RandomNumbersEW4_" + str(Number) + ".txt"
ParamsFile = "/u/flashscratch/j/jmooney3/SimulationsAndParamsABC/OutputParams" + str(Number) + ".txt"

#ABCOutputFile = "/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/ABCResults/ABCOutput_NoError_EW4_" + str(Number) + ".txt"
#Generate and write ABC output files for later
#Change the header based on number of bins based on number of parameter estimates and max number of hets
#ABCneFile = open (ABCOutputFile, 'w')
#ABCneFile.write('SimNeOne\tSimNeTwo\tSimTime\tBIN0\tBIN1\tBIN2\tBIN3\tBIN4\tBIN5\tBIN6\tBIN7\tBIN8\tBIN9\tBIN10\tBIN11\tBIN12\tBIN13\tBIN14\tBIN15\tBIN16\tBIN17\tBIN18\tBIN19\tBIN20\tBIN21\tBIN22\tBIN23\tBIN24\tBIN25\n')
#ABCneFile.close()

#Read through file and skip the header line
f = open (LengthFile, 'r')
for line in f:
	if Reps > 0:
		Split = line.split("\t")
		Bp = int(Split[3])	
		CalledSites.append(Bp)
	Reps = Reps + 1

Reps = Reps - 1

#Open Files
RandFile = open (RandomNumbersFile, 'w')
ParsFile = open (ParamsFile, 'w')

###All of the priors come from PSMC looking at INDV EW10 which had largest ranges 
#sample uniform from 900 to 10000
NeCurr = random.random()* 9100; #10K-900 
NeCurr = int( NeCurr + 900 + 0.5); 

#sample uniform from 11k to 40K
NeTwo = random.random()* 29000; ##40K-11k
NeTwo = int( NeTwo + 11000 + 0.5);

#sample uniform from 1700 to 2400
Time = random.random()*700 ; #2400-1700
Time = int( Time + 1700 + 0.5 );

StringPars = str(NeCurr) + "\t" + str(NeTwo) + "\t"+ str(Time) + "\n"
ParsFile.write(StringPars)

for i in range(0,len(CalledSites)):
	Theta = 4.0 * float(NeCurr) * MutRate * float(CalledSites[i])
	
	RecRateBp = random.gammavariate(0.1, 0.0000001) ### gamma with parameters alpha = 0.1 and beta = 0.0000001. Mean = 1e-8; Var = 3.162278e-08
	RecRateRegion = 4.0* float(NeCurr)* float(RecRateBp) *  float(CalledSites[i])
	
	GensSecEpoch = float(Time) / (4.0*float(NeCurr)); #get to 2nd Epoch
	FractionSecEpoch = float(NeTwo) / float(NeCurr);

	if (CalledSites[i] > LengthThreshold):
		String = str(Theta) + " " + str(RecRateRegion) + " " + str(CalledSites[i]) + " " + str(GensSecEpoch) + " "  + str(FractionSecEpoch) + "\n"
		RandFile.write(String)

RandFile.close()
ParsFile.close()
