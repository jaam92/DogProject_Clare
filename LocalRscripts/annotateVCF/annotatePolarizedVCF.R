###This script will pull information from VCF file for each individual and pair it with corresponding VEP annotation at each site, then output a new file###
###Input files: VEP annotation file and Polarized VCF with ANCHOM,DERHOM,HET, or MISSING as GT
###Output file: File for each individual and chrom with CHROM, POS, GT, ANNOT, IMPACT, SIFT

#Load Libraries and set working directory
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
setwd("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF")

for (chrom in 1:38) {
  #Load file
  
  annotation = read.delim(file=paste("~/Documents/DogProject_Clare/LocalRscripts/parseVEP/parsedVEP_chr", chrom, "_codingAnnots.txt", sep=""))
  
  vcf = fread(file=paste("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AssignGT_WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr",chrom, ".vcf", sep="")) %>% select(-c(V1)) #remove the null column that get's added
  
  #Change column names
  names(vcf)[1] = "CHROM"
  names(vcf)[2] = "POS"
  colnames(vcf) = gsub(".*]","", colnames(vcf)) #remove the [int] from column names
  
  #Add annotations
  annotatedVCF = vcf %>%
    filter(POS %in% annotation$POS) %>%
    mutate(ANNOT = annotation$ANNOT[match(POS,annotation$POS)],
           IMPACT = annotation$IMPACT[match(POS,annotation$POS)],
           SIFT = annotation$SIFT[match(POS,annotation$POS)])
  #Grab individuals
  indivs = colnames(annotatedVCF[,5:79])
  
  for (i in indivs){
    indivAnnotDF = annotatedVCF %>% 
      select(CHROM, POS, i, ANNOT, IMPACT, SIFT) #select columns of interest
    indivID = gsub(":GT", "", i) #clean up id column
    write.table(indivAnnotDF, file=paste("annotatedVEP_" ,indivID, "_Chr" ,chrom, ".txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  
}
    