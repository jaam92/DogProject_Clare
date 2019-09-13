###This script will pull information from VCF file for each individual and pair it with corresponding VEP annotation at each site, then output a new file###
###Input files: VEP annotation file and Polarized VCF with ANCHOM,DERHOM,HET, or MISSING as GT
###Output file: File for each individual and chrom with CHROM, POS, GT, ANNOT, IMPACT, SIFT

#Load Libraries
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)

setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT")

for (chrom in 1:38){

  #Load file
  annotation = read.table(file=paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/VEP/parseVEP/parsedVEP_chr", chrom, "_codingAnnots.txt", sep=""))
  vcf = fread(file=paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RecodedVCF/AssignGT_WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr", chrom, ".vcf", sep="")) %>% select(-c(V1)) #remove the null column that get's added
  
  #Change column names
  names(vcf)[1] = "CHROM"
  names(vcf)[2] = "POS"
  colnames(vcf) = gsub(".*]","", colnames(vcf))
  
  #Add annotations
  annotatedVCF = vcf %>%
    filter(POS %in% annotation$V2) %>%
    mutate(ANNOT = annotation$V3[match(POS,annotation$V2)],
           IMPACT = annotation$V4[match(POS,annotation$V2)],
           SIFT = annotation$V6[match(POS,annotation$V2)])
  #Grab individuals
  indivs = colnames(annotatedVCF[,5:79])

    for (i in indivs){
      indivAnnotDF = annotatedVCF %>%
                      select(CHROM, POS, i, ANNOT, IMPACT, SIFT)
      indivID = gsub(":GT", "", i)
      write.table(indivAnnotDF, file=paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT/annotatedGTwithVEP_",indivID, "_Chr" ,chrom, ".txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }

}

