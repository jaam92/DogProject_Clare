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
RandomNumbersFile = "/u/flashscratch/j/jmooney3/SimulationsAndParamsABC_3epoch/RandomNumbersEWsamps_" + str(Number) + ".txt"
ParamsFile = "/u/flashscratch/j/jmooney3/SimulationsAndParamsABC_3epoch/OutputParams" + str(Number) + ".txt"

#RandomNumbersFile = "testRanNum" + str(Number) + ".txt"
#ParamsFile = "testParam" + str(Number) + ".txt"


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

#Sample uniform from 10 to 1000
NeCurr = random.random()* 990; #1K-10
NeCurr = int( NeCurr + 10 + 0.5); 

#Sample uniform from 10 to 1000
TimeOne = random.random() * 990 ; #1000-10
TimeOne = int( TimeOne + 10 + 0.5 );

#Sample uniform from 1000 to 60K
NeTwo = random.random()* 50000; #60K-1K
NeTwo = int( NeTwo + 1000 + 0.5);

#Sample uniform from 1500 to 15K
TimeTwo = random.random() * 13500 ; #15K-1500
TimeTwo = int( TimeTwo + 1500 + 0.5 );

#Sample uniform from 40K to 150K
NeThree = random.random()* 110000; #150K-40K
NeThree = int ( NeThree + 40000 + 0.5 );


String = str(NeCurr) + "\t" + str(NeTwo) + "\t" + str(TimeOne) + "\t" + str(NeThree) + "\t" + str(TimeTwo) + "\n"
ParsFile.write(String)

for i in range(0,len(CalledSites)):
	Theta = 4.0 * float(NeCurr) * MutRate * float(CalledSites[i])
	
	RecRateBp = random.gammavariate(0.1, 0.0000001) ### gamma with parameters alpha = 0.1 and beta = 0.0000001. Mean = 1e-8; Var = 3.162278e-08
	RecRateRegion = 4.0* float(NeCurr)* float(RecRateBp) *  float(CalledSites[i])
	
	GensOne = float(TimeOne) / (4.0*float(NeCurr));
	FractionOne = float(NeTwo) / float(NeCurr);

	GensTwo = float(TimeTwo) / (4.0*float(NeCurr));
	FractionTwo = float(NeThree) / float(NeCurr);

	if (CalledSites[i] > LengthThreshold):
		String = str(Theta) + " " + str(RecRateRegion) + " " + str(CalledSites[i]) + " " + str(GensOne) + " "  + str(FractionOne) + " " + str(GensTwo) + " " + str(FractionTwo) + "\n"
		RandFile.write(String)

RandFile.close()
ParsFile.close()
