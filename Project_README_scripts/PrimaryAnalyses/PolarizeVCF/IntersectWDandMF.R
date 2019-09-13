####This script takes two input files and generates a single ouput file
  ##It will iterate through all 38 chroms
####Input file 1 is the WD allele file (you can just use the file that I have)
  ###All I did was pull the 1-based pos and WD column from Clare's WD data
####Input file 2 is the file you want to compare against (so change this path)
  ###Print the position, ref, and alt column from your vcf without the header line
####Output file contains that position and ref allele for useable wild dog sites
  ###change the path for your output file as  well 

####To run script on hoffman
  ###module load R
  ###Rscript IntersectWDandMF.R

#Load Libraries
library(dplyr)
library(data.table)
library(tidyr)


for (CHROM in 1:38){
 
#Input file 1: for wild dog you should be able to use this file since I grabbed it from Clares directory 
#All you have to do add to the beginning of the path
  WD = fread(file=paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/WildDog/WD_chr", CHROM , "_CHROM_POS_allele.txt", sep="")) 
  colnames(WD) = c("POS", "Allele")

#Input file 2: for the file that you want to compare against (so change this path to your input file directory)
#Print the pos, ref, and alt alle from the vcf file and remove the header line
  masterFile = read.table(file=paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/WildDog/CanFam3.1_ClareProj/masterFile_CanFam3.1_chr", CHROM ,"_REFALT.txt", sep="")) 
  colnames(masterFile) = c("POS", "REF", "ALT")
  
  #Grab sites that are in both file sets
  df = WD %>%
    filter(POS %in% masterFile$POS) %>% #filter down to only sites in master file
    separate(Allele, c("A1", "A2"), fill = "right") %>% #split column to check hets
    mutate(mf_REF = masterFile$REF[match(POS, masterFile$POS)],
           mf_ALT = masterFile$ALT[match(POS, masterFile$POS)]) #add the REF and ALT from Master File
  
  #Delete the wild dog file to save memory
  rm(WD)  
  
  #Check if we can keep any of the heterozygous sites
  KeepHetDF = df %>%
    filter(!is.na(A2)) %>% #subset hets
    mutate(checkA1_REF = ifelse(A1 == mf_REF, "1", "0"),
           checkA2_ALT = ifelse(A2 == mf_ALT, "1", "0")) %>% #check whether REF or ALT match the masterfile allele
    filter(checkA2_ALT == 0 & checkA1_REF == 1) %>% #if the ALT does not match but the reference does,then we know the REF is ancestral
    select(POS, A1)
  
  #Filter down to sites we can keep 
  Sites2Keep = df %>% 
    filter(is.na(A2)) %>%
    select(POS, A1) %>% #select the position and REF allele
    rbind(KeepHetDF) %>%
    arrange(POS) #order by position
  
  #Write out the position and the REF allele from the wild dog file
  #You will want to change this path to whatever directory you want the output in
  write.table(Sites2Keep, file=paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/WildDog/WD_Alleles_inMasterFile_Chr", CHROM, ".txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  }
  


