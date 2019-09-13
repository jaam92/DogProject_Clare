# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:43:58 2015

@author: christian

Parses Clares Drosophila snpEff file

Format of input file effect string: 
Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon [ | ERRORS | WARNINGS ] )

first line:
['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'ZI103', 'ZI104', 'ZI10', 'ZI112N',...

GOAL:
Write a python script that also checks if the exon is of protein_coding type, and that also reports all genes of possibly overlapping exons, not just the first one.

Output:
chr pos ref alt AC AN effect annGenes allGenes
2L      39466   C       .       0       194     EXON    Cda5    Cda5
2L      39467   G       .       0       194     EXON    Cda5    Cda5
2L      39468   T       C       1       194     NON_SYNONYMOUS  Cda5    Cda5
2L      39469   C       .       0       194     EXON    Cda5    Cda5
2L      39470   T       .       0       194     EXON    Cda5    Cda5
...

If there are several genes annotated for the same position (e.g. overlapping exons), they will be separated by a colon.

Usage (chromosome 2L):
python parse_drosoph_snpEff_file.py BDGP5_75_INFOadded_MERGED_2L_MASKED_bcftools.vcf BDGP5_75_INFOadded_MERGED_2L_MASKED_bcftools_parseCH.txt


@jazlyn altered this script (11/8/2018) because the EFF= annotation is not in the same column every time thus you have to search for it also the new format is ANN= 

"""

import sys, os
import re
import gzip 

fnIN = sys.argv[1]
fnOUT = sys.argv[2]

if len(sys.argv) !=3:
	print "Usage: python parse.py <input vcf> <annotated output file>"
	exit(1)

fileIN = gzip.open(fnIN, "r")
fileOUT = open(fnOUT, "w")

print >> fileOUT, "chr pos ref alt effect annGenes allGenes"

for line in fileIN:
    if line.startswith("#"):
        pass
    else: 
        line = line.split("\t")
        chro = line[0]
        pos = int(line[1])
        ref = line[3]
        alt = line[4]
        info = line[7].split(";")

        EFFIndex = [i for i, s in enumerate(info) if 'ANN=' in s] #find the EFF column
        for x in EFFIndex:
            num = int(x) #convert the items in the list to integers
            EFF = info[x].split("=")[1].split(",")
            EFF_Effect = [x.split("|")[1] for x in EFF]
            EFF_String = [x.split("|") for x in EFF]

            #EFF_EXONS = [[j, i] for i,j in zip(EFF_String, EFF_Effect) if re.search("^EXON.*",  j) and i[5] == "protein_coding"]  # check if there are exons that are also protein coding (there are some ncRNA exons...)
            EFF_SYN = [[j, i] for i,j in zip(EFF_String, EFF_Effect) if re.search("^synonymous*",  j) and i[7] == "protein_coding"]
            #print(EFF_SYN)  
            EFF_NONSYN = [[j, i] for i,j in zip(EFF_String, EFF_Effect) if re.search("^missense.*",  j) and i[7] == "protein_coding"]
            #print(EFF_NONSYN)  
            #EFF_INTRON = [[j, i] for i,j in zip(EFF_String, EFF_Effect) if re.search("^INTRON.*",  j) and i[5] == "protein_coding"]  
        
            
            if EFF_SYN == [] and EFF_NONSYN == []:
                continue  # the site is not in an exon
                
            EFF_GENE = []
            if EFF_NONSYN != []: EFF_GENE.append(EFF_NONSYN)
            if EFF_SYN != []: EFF_GENE.append(EFF_SYN)
            #if EFF_EXONS != []: EFF_GENE.append(EFF_EXONS)   # Note that exons technically contain UTR's, but in this file they do not! I.e. an exon position is never an UTR, it is always coding.
            allGenes = set([x[0][1][4] for x in EFF_GENE])   # all annotated genes at this position
            
            if EFF_NONSYN != []:
                effect = "NON_SYNONYMOUS"
                annGenes = set([x[1][4] for x in EFF_NONSYN])   # genes that lead to the nonsynonymous SNP
            elif EFF_SYN != []:
                effect = "SYNONYMOUS"
                annGenes = set([x[1][4] for x in EFF_SYN])   # genes that lead to the synonymous SNP
            #elif EFF_EXONS != []:
            #    effect = "EXON"
            #    annGenes = set([x[1][4] for x in EFF_EXONS])   # genes that lead to the exon
              
            # print EFF_SYN
            # print EFF_EXONS
            # print EFF_SYN
            print "\t".join([chro, str(pos), ref, alt, effect, ",".join(annGenes), ",".join(allGenes)])
            print >> fileOUT, "\t".join([chro, str(pos), ref, alt, effect, ",".join(annGenes), ",".join(allGenes)])
