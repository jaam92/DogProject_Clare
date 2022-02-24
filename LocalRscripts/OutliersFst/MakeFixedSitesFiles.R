#Load Libraries
library(data.table)
library(tidyverse)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/")

indivsEW = read.table(file = "~/Documents/DogProject_Clare/LocalRscripts/IndivFiles/EW_indivList_downsample_N9.txt")
indivsTM = read.table(file = "~/Documents/DogProject_Clare/LocalRscripts/IndivFiles/TM_indivList_downsample_N9.txt")
indivs = rbind(indivsEW, indivsTM)


#Empty data frame to fill with summary data
summaryInfo = data.frame()

for (chrom in 1:38){
  vcf = fread(file=paste("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AssignGT_WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr", chrom, ".vcf", sep="")) %>% select(-c(V1)) #remove the null column that get's added
  
  #Change column names
  colnames(vcf) = gsub(".*]","", colnames(vcf))
  colnames(vcf) = gsub(":GT", "", colnames(vcf))
  N6_indivs = colnames(vcf)[which(colnames(vcf) %in% indivs$V1)]
  
  GTCounts = vcf %>%
    select(CHROM, POS, N6_indivs) %>%
    filter_all(all_vars(. != "Missing")) %>%
    filter_all(all_vars(. != "Het"))
  
  rm(vcf)
  
  #Replace pattern with value
  GTCounts[GTCounts=="DerHom"]<- 2
  GTCounts[GTCounts=="AncHom"]<- 0
  
  #Count the number of DERHOM or ANCHOM 
  EW = GTCounts %>%
    select(starts_with("EW")) %>%
    mutate_all(funs(parse_number(.))) %>% 
    rowSums()
  
  #Count the number of DERHOM or ANCHOM 
  TM = GTCounts %>%
    select(starts_with("TM")) %>% 
    mutate_all(funs(parse_number(.))) %>%
    rowSums()
  
  #keep those sites that are fixed as opposites in TM and EW
  combo = GTCounts %>%
    select(CHROM, POS) %>%
    cbind.data.frame(EW,TM) %>%
    filter(EW == "0" & TM == "18" | EW == "18" & TM == "0")
  
  
  #Data frame with info from all chromosomes
  summaryInfo = bind_rows(summaryInfo, combo)

}

write.table(summaryInfo, file = "FixedSites_EWvsTM_N9.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


