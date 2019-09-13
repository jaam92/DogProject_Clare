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
#sample uniform from 10 to 10000
NeCurr = random.random()* 9990; #10K-10 
NeCurr = int( NeCurr + 10 + 0.5); 

#sample uniform from 11k to 70K
NeTwo = random.random()* 59000; ##70K-11k
NeTwo = int( NeTwo + 11000 + 0.5);

#sample uniform from 10 to 15000
Time = random.random()*14990 ; #15000-10
Time = int( Time + 10 + 0.5 );

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
