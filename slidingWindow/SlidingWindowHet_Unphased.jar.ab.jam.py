# script to count number of genotypes called and number of heterozygotes per sample
# input file is a VCF file that has been filtered
# supply input file, window size, step size, chromosome number
# example: python ./SlidingWindowHet.py input.vcf.gz 100000 10000 1

import sys
import pysam
import os
import gzip
import argparse

def parse_args():
	"""
	Parse command-line arguments
	"""
	parser = argparse.ArgumentParser(description="This script computes sliding window heterozygosity.")

	parser.add_argument(
            "--vcf", required=True,
            help="REQUIRED. input vcf file, this must have passed through Step 11 so *NO* bad sites remain in the file")
	parser.add_argument(
			"--window_size", required=True,
			help="REQUIRED. window size ")
	parser.add_argument("--step_size", required=True,
			help="REQUIRED. Name of output file.")
	parser.add_argument("--chromNum", required=True,
			help="REQUIRED. chromosome of query VCF")
	args = parser.parse_args()
	return args
args=parse_args()

# open input file (gzipped VCF file), make sure the VCF file is indexed (if not, create index)
filename = args.vcf
VCF = gzip.open(filename, 'r')

if not os.path.exists("%s.tbi" % filename):
	pysam.tabix_index(filename, preset="vcf")
parsevcf = pysam.Tabixfile(filename)


# set variables
window_size = int(args.window_size)
step_size = int(args.step_size)
chromo = args.chromNum
het = ['0/1','1/0'] 

#Chromosome size humans
#chromo_size={'1':248956422,'2':242193529,'3':198295559,'4':190214555,'5':181538259,'6':170805979,'7':159345973,'8':145138636,'9':138394717,'10':133797422,'11':135086622,'12':133275309,'13':114364328,'14':107043718,'15':101991189,'16':90338345,'17':83257441,'18':80373285,'19':58617616,'20':64444167,'21':46709983,'22':50818468}
#Chromosome size dogs
chromo_size={'chr1':122678785,'chr2':85426708,'chr3':91889043,'chr4':88276631,'chr5':88915250,'chr6':77573801,'chr7':80974532,'chr8':74330416,'chr9':61074082,'chr10':69331447,'chr11':74389097,'chr12':72498081,'chr13':63241923,'chr14':60966679,'chr15':64190966,'chr16':59632846,'chr17':64289059,'chr18':55844845,'chr19':53741614,'chr20':58134056,'chr21':50858623,'chr22':61439934,'chr23':52294480,'chr24':47698779,'chr25':51628933,'chr26':38964690,'chr27':45876710,'chr28':41182112,'chr29':41845238,'chr30':40214260,'chr31':39895921,'chr32':38810281,'chr33':31377067,'chr34':42124431,'chr35':26524999,'chr36':30810995,'chr37':30902991,'chr38':23914537,'chrX':123869142}

# get list of samples
samples=[]
for line in VCF:
	if line.startswith('##'):
		pass
	else:
		for i in line.split()[9:]: samples.append(i)
		break


# get first and last positions in chromosome
for line in VCF:
	if line[0] != '#':
		start_pos = int(line.strip().split()[1])
		end_pos = chromo_size[chromo]
		break


# create output file
output = open(filename + '_het_%swin_%sstep.txt' % (window_size, step_size), 'w')
output.write('chromo\twindow_start\tsites_total\tsites_passing\tcalls_%s\tmissing_%s\thets_%s\thomRef_%s\thomAlt_%s\n' % ('\tcalls_'.join(samples), '\tmissing_'.join(samples), '\thets_'.join(samples), '\thomRef_'.join(samples), '\thomAlt_'.join(samples)) )


# Fetch a region, ignore sites that fail filters, tally total calls and heterozygotes		
def snp_cal(chromo,window_start,window_end):

	rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (chromo, window_start, window_end), parser=pysam.asTuple()))
	
	sites_total,sites_passing=0,0
	calls=[0]*len(samples)
	hets=[0]*len(samples)
	missing=[0]*len(samples)
	homRef=[0]*len(samples)
	homAlt=[0]*len(samples)

	for line in rows:
		sites_total+=1	
		sites_passing+=1
		for i in range(0,len(samples)):
			GT=line[i+9]	
			if GT[:1]!='.': calls[i]+=1	
			if (GT[:3] in het): hets[i]+=1
			elif GT[:3]=='./.': missing[i]+=1
			elif GT[:3]=='0/0': homRef[i]+=1
			elif GT[:3]=='1/1': homAlt[i]+=1
			else: print GT[:3]

	output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (chromo,window_start,sites_total,sites_passing,'\t'.join(map(str,calls)),'\t'.join(map(str,missing)),'\t'.join(map(str,hets)),'\t'.join(map(str,homRef)),'\t'.join(map(str,homAlt))) )

	
# initialize window start and end coordinates
window_start = start_pos
window_end = start_pos+window_size-1


# calculate stats for window, update window start and end positions, repeat to end of chromosome
while window_end <= end_pos:	
			
	if window_end < end_pos:
		snp_cal(chromo,window_start,window_end)

		window_start = window_start + step_size
		window_end = window_start + window_size - 1

	else:
		snp_cal(chromo,window_start,window_end)
		break	
		
else:
	window_end = end_pos
	snp_cal(chromo,window_start,window_end)


# close files and exit
VCF.close()
output.close()

exit()
