import os
import sys

if len(sys.argv) != 4:
	print "Usage: python clean_vcf.py <WildDog_pos_allele file> <vcf file> <output file>"
	exit(1)

WildDogDataFilename = sys.argv[1]
seqDataFilename = sys.argv[2]
recodedSeqDataFilename = sys.argv[3]

goodValues = ['A','T','C','G'] #check to see that REF and ALT positions are valid alleles

#only do work if both files exist
if (os.path.isfile(WildDogDataFilename) and os.path.isfile(seqDataFilename)):
    
	#create dictionary with rs id and allele
	with open(WildDogDataFilename, 'r') as f:
		WildDogData = dict(line.strip().split(None, 1) for line in f)
 
	#open vcf file 
	with open(seqDataFilename, 'r') as datafile:
		with open(recodedSeqDataFilename, 'w') as outFile:
			for line in datafile:
				line = line.lstrip()
				#skip lines that start with # character  
				if not line.startswith("#"):
					line = line.split('\t')
					#print line
					position = line[1]
					CanidRef = line[2]                        
					CanidAlt = line[3]                        				    
					recodedLine = ''

					#only work with lines if the position is in WildDog data 
					if position in WildDogData:
						k = WildDogData[position] 
						#print k
						if (CanidAlt in goodValues and CanidRef in goodValues and k in goodValues):	#only process lines with A,T,C,G as value
							if CanidRef == k:
								outFile.writelines("\t".join(line))
								#print 'REF matched'
				
							elif CanidAlt == k:			
								#print 'ALT matched'
								recodedLine = ''
								
								#write initial columns that don't hold people data but stop at ref value
								for dataCol in line[:2]:
									if recodedLine != '':
										recodedLine = recodedLine + '\t' + dataCol#each col tab separated
									else:
										recodedLine = dataCol#starting col by itself
			    
								#write the new reference value and alternate value from dictionary
								recodedLine = recodedLine + '\t' + k + '\t' + CanidRef

								#write remaining initial columns that don't hold people data
								#for dataCol in line[5:9]:
								#	recodedLine = recodedLine + '\t' + dataCol#each col tab separted

								#write the people data starting at column idx 9 & check length of genotype col
								for dataCol in line[4:]:
									if len(dataCol) > 4:
										print '** ERROR genotype col is > 4'
										print dataCol
										sys.exit('*** ERROR exiting ***')
									# strip if extra end of line
									if len(dataCol) == 4:
										#print 'len', len(dataCol)
										#print '['+dataCol+']'
										dataCol=dataCol.strip()
										# check fix worked
										if len(dataCol) !=3:
											sys.exit('** EXITING stripping of end of line didnt work there is still error **')

								#start recoding columns
									if dataCol == '1/1':
										newCol = '0/0'
									elif dataCol == '1/0':
										newCol = '0/1'
									elif dataCol == '0/1':
										newCol = '1/0'
									elif dataCol == '0/0':
										newCol = '1/1'
									elif dataCol == './.':
										newCol = './.'
									else:
										print dataCol
										sys.exit('**exiting error**')
										
									recodedLine = recodedLine + '\t' + newCol #each col tab separted and add end of line before writing
								
								recodedLine = recodedLine + '\n' #This will add line return after entire recode line is complete
								
								outFile.writelines(recodedLine)#write newly constructed line
								
							else:
								print 'no match', position, k, CanidRef, CanidAlt #if WildDog Ref doesn't match either REF or ALT skip processing 
									
						else:
							print 'ERROR INVALID ALLELE at REF or ALT or WildDog', position, CanidRef, CanidAlt, k #if allele isn't A,T,C,G skip processing 
							#outFile.writelines(line) #keep line even if there was error in allele value 
				else:
					outFile.writelines(line) #keeping initial hashed lines



			
