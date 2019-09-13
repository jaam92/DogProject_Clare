###This script will identify positions fall within a given range by appending a column to the input file with either 1 or 0
###The input file needs to contain column with positions and first line should be a header line
###The first column in the range file should be the start and the second column should be the end 
###The output file will contain the original information from the input file and a new column named after the given range descriptor 

import sys
import os
import itertools

numberOfArgs = len(sys.argv)

if (numberOfArgs < 7):
	print "Missing arguments.   Usage: python addPositionAttributes.py inputFilename posColumnNumber outputFilename [ranges|positions|positionAnnotation] range1Filename range1Descriptor"
	exit(1)
if (numberOfArgs > 7 and ((numberOfArgs-1)%2)!=0):
	print "Missing arguments.   Usage: python addPositionAttributes.py inputFilename posColumnNumber outputFilename [ranges|positions|positionAnnotation] range1Filename range1Descriptor range2Filename range2Descriptor"
	exit(1)

numberOfRangeFiles = (numberOfArgs - 5)/2

inputFilename = sys.argv[1]
posColumnIndex = int(sys.argv[2]) - 1
outputFilename = sys.argv[3]
typeOfComparison = sys.argv[4]
intermediateFile = "chr15to38outResults.tmp"
endOfLine = '\n'
tab = '\t'

if (typeOfComparison != "ranges" and typeOfComparison != "positions" and typeOfComparison != "positionAnnotation"):
	print "Unknown comparison type.  Type must be specified as ranges, positions or positionAnnotation.  Example: python addPositionAttributes.py inputFilename posColumnNumber outputFilename ranges range1Filename range1Descriptor"
	exit(1)

##functions to be used in main loop
# getPostion function has 2 parameters:
#	dataline is the whole line from dat file
#	positionColumn is the zero-based column index of the position values in the dat file
def getPosition(dataline, positionColumn):
	dataColumns = dataline.strip().split(tab)
	pos = dataColumns[positionColumn]
	return int(pos)

def getPositionAnnotationDictionary(positionsFilename):
	positionDict = {}
	if (os.path.isfile(positionsFilename)):
		with open(positionsFilename, 'r') as inputFile:
			positionDict = dict(line.strip().split(None, 1) for line in inputFile);
	return positionDict

	
def getPositionsList(positionsFilename):
	positions = []
	if (os.path.isfile(positionsFilename)):
		with open(positionsFilename, 'r') as inputFile:
			allLines = inputFile.read();
			try:
				positions = [int(x) for x in allLines.split()]
			except:
				pass  #do nothing; values won't be included in annotation positions
	return positions

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

	
#main loop
for i in range(1, numberOfRangeFiles+1):
	rangeFilename = sys.argv[2*i+3]
	rangeColumnName = sys.argv[2*i+4]
	intermediateFilename = str(i) + "_" + intermediateFile
	headerLineIsDone = False
	
	possiblePositions = [] #initializer
	if (typeOfComparison == "ranges"):
		posRanges = getRanges(rangeFilename)
		#print posRanges
		from itertools import chain
		possiblePositions = set(chain(*(xrange(start, end+1) for start, end in posRanges)))
	if (typeOfComparison == "positions"):
		possiblePositions = getPositionsList(rangeFilename)
	if (typeOfComparison == "positionAnnotation"):
		possiblePositionsDict = getPositionAnnotationDictionary(rangeFilename)
	#print possiblePositions
	#print possiblePositionsDict
	
	#if this is first pass then use specified input filename
	if (i == 1):
		loopInputFilename = inputFilename
	else:
		loopInputFilename = loopOutputFilename
	
	#assume multiple range files and so output is defaulted as intermediate filename
	loopOutputFilename = intermediateFilename
	#unless we are on the last run, then use specified output filename
	if (i == numberOfRangeFiles):
		loopOutputFilename = outputFilename
		
	#read thru input data and write out attribute flag for each position
	with open(loopOutputFilename, 'w') as outputFile:
		with open(loopInputFilename, 'r') as inputFile:
			for line in inputFile: #read next input line
				if headerLineIsDone:
					posValue = getPosition(line, posColumnIndex) #since position is in the 1st column then send 0 as the index
					#print posValue
					
					if (typeOfComparison != "positionAnnotation"):
						valueExists = "0"
						if posValue in possiblePositions:
							valueExists = "1"
					else:
						if possiblePositionsDict.has_key(str(posValue)):
							valueExists = possiblePositionsDict[str(posValue)]
						else:
							valueExists = "NA"
							
					outputFile.write(line.rstrip() + tab + valueExists + endOfLine)
				else:
					outputFile.write(line.rstrip() + tab + rangeColumnName + endOfLine)
					headerLineIsDone = True
	
	#clean up memory
	if (typeOfComparison == "positionAnnotation"):
		del possiblePositionsDict
	
	possiblePositions = []
	#clean up disk space
	if (i != 1): #only delete intermediate files after the 1st loop (original input)	
		os.remove(loopInputFilename)

print "Processing Completed. Output File: ", outputFilename
