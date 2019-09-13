import sys 
import random 
import math 
import numpy as np

#Input args
#for arg in sys.argv:
	#print arg
EmpiricalFile = sys.argv[1] #input file for empirical data
LengthThreshold = int(sys.argv[2]) #Threshold cutoff for simulation 200bp
SimData = sys.argv[3] #input file for simulated data
TaskID = sys.argv[4] #SGE task id
BinomialFlag = int(sys.argv[5]) 
ErrorRate = float(sys.argv[6])
CallableSitesEmpirical = int(sys.argv[7]) #col number with the number of Callable sites
BinsToTake = int(sys.argv[8]) #minimum bins
UpperBinLimit = int(sys.argv[9]) #maximum bins (M)
ParamFileDirectory = sys.argv[10] #provide the directory that contains parameter files
OutputFileHeader = sys.argv[11] #provide the directory and file header for ABC output

#Output files
ParamsFile = "OutputParams" + TaskID + ".txt"
ErrorFile = "OutputParamsError" + TaskID + ".txt"
ABCOutputFile = OutputFileHeader + TaskID + ".txt"

#Set up Counters and Lists
Reps = 0
ErrorRate = 0
SimTaskID = 0
AbsoluteDistance = 0
AbsoluteDistanceRob2016 = 0
HeterozygotesEmipircal = []
NumberOfBases = []
SimHeterozygotes = []
SimPi = []
CountsSim = []
CountsData = []

#Iterate through emipircal data
f = open (EmpiricalFile,'r')
for line in f:
	TLine = line.rstrip('\n')
	if ( Reps > 0 ): #skip header line
		Split = TLine.split("\t")
		CallableSites = int(Split[CallableSitesEmpirical])
		HetCol = Split[CallableSitesEmpirical+8]	
		if (CallableSites > LengthThreshold):
			HeterozygotesEmipircal.append(int(HetCol))
			NumberOfBases.append(CallableSites)
	Reps = Reps + 1
Reps = Reps - 1
f.close()

#Iterate through simulated data
k = open (SimData,'r')
if BinomialFlag == 1:
	ParamsFile = ParamFileDirectory + str(ErrorFile)
else:
	ParamsFile = ParamFileDirectory + str(ParamsFile)	

for line in k:
	TLine = line.rstrip('\n')
	Split = TLine.split("\t")
	PiColSim = Split[0]
	HetColSim = Split[1]
	#print PiColSim,HetColSim
	if BinomialFlag == 1:
		Errors = int(np.random.binomial(NumberOfBases[SimTaskID],ErrorRate,1))
		SimHeterozygotes.append(int(HetColSim) + Errors)
	else:
		SimHeterozygotes.append(int(HetColSim))
		SimPi.append(float(PiColSim))	
	SimTaskID = SimTaskID + 1
k.close()

#Number of lines in each empirical data and simulations (should be equal)
LengthOne = len(HeterozygotesEmipircal)
LengthTwo = len(SimHeterozygotes)
#print "Number Sims = " + str(LengthOne) + " " + str(LengthTwo)

#Create list to store heterozygous bins
for i in range(0, 4000): #Maximum amount of hets is 4000 since windows are 1Kb and we are using 4 individuals
	CountsSim.append(0)
	CountsData.append(0)

#Make a list of lists with counts of hets per x bin
for i in range(0, LengthTwo): 
	CountsSim[SimHeterozygotes[i]] = CountsSim[SimHeterozygotes[i]] + 1
	CountsData[HeterozygotesEmipircal[i]] = CountsData[HeterozygotesEmipircal[i]] + 1
#print CountsSim

#Calculate Distance metric from Robinson 2016 
for i in range(0,BinsToTake):
	#AbsoluteDistance = AbsoluteDistance + float(i) * float(abs(float(CountsSim[i]) - float(CountsData[i])))
	#AbsoluteDistanceRob2016 = AbsoluteDistanceRob2016 + float(abs(float(CountsSim[i]) - float(CountsData[i])))
#uncomment to calculate the distance and normalized by the empirical data
	if CountsData[i] == 0:
		AbsoluteDistance = AbsoluteDistance + float(i) * 0
		AbsoluteDistanceRob2016 = AbsoluteDistanceRob2016 + 0
	else:
		AbsoluteDistance = AbsoluteDistance + float(i) * float(abs(float(CountsSim[i]) - float(CountsData[i]))/(float(CountsData[i])))
		AbsoluteDistanceRob2016 = AbsoluteDistanceRob2016 + float(abs(float(CountsSim[i]) - float(CountsData[i]))/(float(CountsData[i])))
	#print "Distance in " + str(i) + " = " + str(CountsSim[i]) + " and " + str(CountsData[i]) 
		print AbsoluteDistanceRob2016
#Write sim parameters, distance metrics, and sim hets per bin starting at 0 going to UpperBin to output
ParFile = open (ParamsFile , 'r')
ABCFile = open (ABCOutputFile, 'a')
for line in ParFile:
	TLine = line.rstrip('\n')
	ABCFile.write(TLine + "\t" + str(AbsoluteDistance) + "\t" + str(AbsoluteDistanceRob2016))
	for i in range(0,UpperBinLimit+1):
		ABCFile.write("\t" + str(CountsSim[i]))	
	ABCFile.write("\n")

ParFile.close()
ABCFile.close()

