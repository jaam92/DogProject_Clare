#Load Libraries
library(data.table)
library(tidyverse)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/SFS")
pops = c("BC","LB","PG","TM","AW","EW","IR")

for (i in pops){
  indivs = read.table(file = paste("~/Documents/DogProject_Clare/LocalRscripts/IndivFiles/" ,i, "_indivList_downsample_N6.txt", sep = ""))
  
  #Empty data frame to fill with summary data
  summaryInfo = data.frame()
  
  for (chrom in 1:38){
    vcf = fread(file=paste("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AssignGT_WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr", chrom, ".vcf", sep="")) %>% select(-c(V1)) #remove the null column that get's added
    
    #Change column names
    colnames(vcf) = gsub(".*]","", colnames(vcf))
    colnames(vcf) = gsub(":GT", "", colnames(vcf))
    N6_indivs = colnames(vcf)[which(colnames(vcf) %in% indivs$V1)]
    
    #Reformat 
    GTCounts = vcf %>%
      select(N6_indivs) %>%
      filter_all(all_vars(. != "Missing")) %>% #remove any sites with missing data
      group_by_all() %>% 
      count() %>% #count of all patterns in rows
      ungroup()
    
    rm(vcf)
    
    #Replace pattern with value
    GTCounts[GTCounts=="DerHom"]<- 2
    GTCounts[GTCounts=="Het"]<- 1
    GTCounts[GTCounts=="AncHom"]<- 0
    names(GTCounts)[7] = "Count"
    
    #count of the frequency of each patter
    FreqBin = GTCounts %>% 
      select(-c(Count)) %>% 
      mutate_all(funs(parse_number(.))) %>% 
      rowSums() 
    
    #Data frame with frequency and counts
    Polarized = cbind.data.frame(FreqBin, GTCounts$Count) %>% 
      group_by(FreqBin) %>% 
      summarise(sum = sum(`GTCounts$Count`)) %>%
      mutate(chromosome = chrom)
      
    #Data frame with all chromosomes
    summaryInfo = bind_rows(summaryInfo, Polarized)
  
  }

  write.table(summaryInfo, file = paste(i, "_SFS_allChroms_SummaryFile_N6.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


}



