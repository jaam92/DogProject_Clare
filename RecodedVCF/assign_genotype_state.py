import sys
import os 

if len(sys.argv) !=3:
	print "Useage: python assign_genotype_state.py <vcf file> <out file>"
	exit(1)

seqDataFilename = sys.argv[1]
recodeGenoStateFilename = sys.argv[2]

if (os.path.isfile(seqDataFilename)):
	with open(seqDataFilename,'r') as inputFile:
		with open (recodeGenoStateFilename,'w') as outFile:
			for line in inputFile:
				line = line.lstrip()
				#skip lines with hash
				if not line.startswith('#'):
					line = line.split('\t')
					recodedLine = ''

					#write out inital lines that are not being recoded
					for dataCol in line[0:4]:
						recodedLine = recodedLine + '\t' + dataCol

					#Check length of genotypes
					for dataCol in line[4:]:
						if len(dataCol) > 4:
							print '** ERROR genotype col is > 4'
							sys.exit('*** ERROR exiting ***')
						# strip if extra end of line
						if len(dataCol) == 4:
							#print 'len', len(dataCol)
							#print '['+dataCol+']'
							dataCol=dataCol.strip()
							# check fix worked
							if len(dataCol) !=3:
								sys.exit('** EXITING still formatting error **')
						
						#start converting to genotype state
						if dataCol == '1/1':
							newCol = 'DerHom'
						elif dataCol == '1/0':
							newCol = 'Het'
						elif dataCol == '0/1':
							newCol = 'Het'
						elif dataCol == '0/0':
							newCol = 'AncHom'
						elif dataCol == './.':
							newCol = 'Missing'
						else:
							print dataCol
							sys.exit('**Exiting Genotype Error**')
										
						recodedLine = recodedLine + '\t' + newCol #each col tab separted and add end of line before writing

					recodedLine = recodedLine + '\n' #This will add line return after entire recode line is complete
					
					outFile.writelines(recodedLine)#write newly constructed line						
				
				else:
					outFile.writelines(line) #write hashed info 
