###This script takes a file that contains ranges only (start col 1 and end col2) and  
###A file that contains positions to be searched. You need to specify what column the positions to be searched are in
###Script will calculate the length of the range, how many positions fall within range, and the proportion of the range that is covered by given sites

import os
import sys

numArgs = len(sys.argv)

if (numArgs < 5):
	print "Missing arguments. Usage: python ROHQualCounts.py exonRangeFile posFile posColumNumber outputFile"
	exit(1)

#Define arguments
RangeFilename = sys.argv[1]
positionsFilename = sys.argv[2]
posColumnIndex = int(sys.argv[3]) - 1
outputFilename = sys.argv[4]
endOfLine = '\n'
tab = '\t'

#Set Counters
COUNT_RANGES=0

#Define Functions
#create a list of ranges from input range file
def getRanges(rangesFilename):
    posRanges = []
    if (os.path.isfile(rangesFilename)):
        with open(rangesFilename, 'r') as inputFile:
            for line in inputFile:
                try:
                    a, b = (int(x) for x in line.split())
                    posRanges.append((a,b))
                except:
                    pass  #do nothing; values won't be included in position ranges
    return posRanges

#getPostion function has 2 parameters to create a list of searchable positions:
#position file is the input position file
#positionColumn is the zero-based column index of the position values in the dat file
def getPosition(posFilename, positionColumn):
	positions = []
	if (os.path.isfile(posFilename)):
		with open(posFilename, 'r') as inputFile:
			for dataline in inputFile:
				dataColumns = dataline.strip().split(tab)
				pos = int(dataColumns[positionColumn])
				positions.append(pos)
	return positions

#Implement Functions
GoodPos = getPosition(positionsFilename, posColumnIndex)
ROHRanges = getRanges(RangeFilename)
print "Gathered information from input files"

#Start Script
print "Starting to evaluate ROH ranges..."
with open(outputFilename, 'w') as outputFile:
	for (lowerBound, upperBound) in ROHRanges:
		LengthROH = upperBound-lowerBound
		possiblePositions = set(range(lowerBound, upperBound+1))#create a set of possible positions within the range
		COUNT_RANGES+=1
		COUNT_POS_IN_RANGE=0
		
		#Check whether any of our searchable positions match the positions within our ROH Range
		for a in GoodPos:
			if a in possiblePositions:
				COUNT_POS_IN_RANGE+=1

			else:
				pass

		possiblePositions = [] #delete the set of positions within range now that we are done using it

		if (int(COUNT_POS_IN_RANGE)!= 0 and int(LengthROH)!=0):
			propGoodSitesinROH = int(COUNT_POS_IN_RANGE)/float(LengthROH)
		else:
			propGoodSitesinROH = 0
		
		outputFile.write(str(lowerBound) + tab + str(upperBound) + tab + str(COUNT_POS_IN_RANGE) + tab + str(LengthROH) + tab + str.format("{0:.3f}", propGoodSitesinROH) + endOfLine) 
		print " Finished with range", COUNT_RANGES   
