####This script will make the permutation annotation files from a single individual's annotation file (you can use any individual bc all VEP and GERP annotations are the same and everyone has the same number of sites annotated) 

#Load files
library(dplyr)

#Read file in 
annots = read.delim(file = "~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AllChroms/annotatedGTwithVEP_AW10_allChroms.txt")
chromLengths = read.delim(file = "~/Documents/DogProject_Clare/LocalRscripts/chromosomeLengths.txt", sep = " ")

#Grab columns I want for regular snps
VEPFile = annots %>%
  select(CHROM, POS, ANNOT) 

#Grab columns I want for snps with GERP scores 
GERPScoreFile = annots %>%
  filter(GERPScore >= 4 | GERPScore <= 2) %>% #Filter to neutral or deleterious GERP Scores 
  select(CHROM, POS, ANNOT,GERPScore) %>%
  filter(GERPScore <= 2 | GERPScore >= 4 & ANNOT != "SY") %>% #Remove SY variation GERP says is deleterious
  mutate(GERPScoreANNOT = ifelse(GERPScore >= 4, "DEL", "NEU"))

#Grab counts of SYN ANNOT that GERP says are Deleterious
CountSYasPutDel = annots %>%
  filter(GERPScore >= 4 & ANNOT == "SY") %>%
  count(CHROM) %>%
  mutate(chromLength = chromLengths$LENGTH[match(CHROM,chromLengths$CHROM)],
         propIncongruentCalls = n/chromLength)

#Write to output files 
write.table(VEPFile, file= "~/Documents/DogProject_Clare/LocalRscripts/ROHPermutation/AnnotationFiles/allChromsVEPAnnotation.txt", row.names = FALSE , col.names= TRUE, quote = FALSE, sep = "\t" )

write.table(GERPScoreFile, file= "~/Documents/DogProject_Clare/LocalRscripts/ROHPermutation/AnnotationFiles/allChromsGERPScoreAnnotation.txt", row.names = FALSE , col.names= TRUE, quote = FALSE, sep = "\t" )

write.table(CountSYasPutDel, file= "~/Documents/DogProject_Clare/LocalRscripts/ROHPermutation/AnnotationFiles/allChromsIncongruentAnnotation.txt", row.names = FALSE , col.names= TRUE, quote = FALSE, sep = "\t" )




