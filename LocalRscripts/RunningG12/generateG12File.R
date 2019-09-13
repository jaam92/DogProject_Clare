#Load Libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyverse)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/RunningG12/")
pops = c("LB","PG","TM","AW","EW")

for (i in pops){
  
  for (chrom in 37:38){
    vcf = fread(file=paste("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AssignGT_WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr", chrom, ".vcf", sep="")) %>% select(-c(V1)) #remove the null column that get's added
    
    #Change column names
    colnames(vcf) = gsub(".*]","", colnames(vcf))
    colnames(vcf) = gsub(":GT", "", colnames(vcf))
    
    #Subset populations of interest and recod
    recodedVCF = vcf %>% 
      mutate_at(vars(starts_with(i)), funs(case_when(. == "AncHom" ~ REF,
                                                     . == "DerHom" ~ ALT,
                                                     . == "Missing" ~ "N",
                                                     . == "Het" ~ "."))) %>%
      select(POS,contains(i)) 
    
    #Remove fixed sites
    end = dim(recodedVCF)[2] 
    numSamps = end - 1 
    keep = apply(recodedVCF[2:end], 1, function(x) length(unique(x[!is.na(x)])) > 1)
    FinalDF = recodedVCF[keep,] 
    
    #Uncomment to remove sites based on amount of missing data (max missing data is 20% because of Clare's filtering)
    #intermidiateDF = recodedVCF[keep,] 
    #intermidiateDF$NumMissing = rowSums(intermidiateDF == "N")
    #FinalDF = intermidiateDF %>%
      #filter(NumMissing < 2) %>%
      #select(-c(NumMissing))
    
    write.table(FinalDF, file = paste(i,"_n", numSamps,"_Chr", chrom ,"_inputG12haps.txt", sep = ""), sep = "," , quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  
}
    
    
    