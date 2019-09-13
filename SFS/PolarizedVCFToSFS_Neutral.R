#Load Libraries
library(data.table)
library(dplyr)
library(stringr)
library(tidyverse)
library(GenomicRanges)

#Load files
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/SFS")

#Load Neutral regions and make 1 based
Neutral = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/CanFam3.1_ExonIntronNeutral/NeutralRegions/allChroms_outside.ensemble_rm.phastConsEl_cutoff_0.4.bed") %>%
  mutate(START1based = START+1)

#Iterate through populations
pops = c("BC","LB","PG","TM","AW","EW","IR")
for (i in pops){
  indivs = read.table(file = paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/" ,i, "_indivList_downsample_N6.txt", sep = ""))
  
  #Empty data frame to fill with summary data
  summaryInfo = data.frame()
  
  for (chrom in 1:38){
    vcf = fread(file=paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RecodedVCF/AssignGT_WildDogasREF_reformatted_masterFile_allCanidsN75_ANmt135_Chr", chrom, ".vcf", sep="")) %>% select(-c(V1)) #remove the null column that get's added

#Change column names
    colnames(vcf) = gsub(".*]","", colnames(vcf))
    colnames(vcf) = gsub(":GT", "", colnames(vcf))
    vcf$POS_END = vcf$POS
    N6_indivs = colnames(vcf)[which(colnames(vcf) %in% indivs$V1)]
    
    #Pull neutral regions for specific chrom
    chromo = paste("chr",chrom, sep = "")
    subsetNeutral = Neutral %>%
      filter(CHROM == chromo) %>%
      select(CHROM, START1based, END)
    
    #Make GRanges object to find overlaps
    posRanges = with(vcf, GRanges(CHROM, IRanges(start=POS, end = POS_END)))
    neutralRanges = with(subsetNeutral, GRanges(CHROM, IRanges(start=START1based, end = END)))
    
    #Find Positions Overlapping Neutral
    POSNeutralOverlaps = findOverlaps(query = posRanges, subject = neutralRanges, type = "within") %>%
      as.data.frame() %>%
      distinct(queryHits)
    
    #Keep Positions within Neutral
    withinNeutral = vcf[c(POSNeutralOverlaps$queryHits),]
	
    rm(vcf)
    
    #Reformat 
    GTCounts = withinNeutral %>%
      select(N6_indivs) %>%
      filter_all(all_vars(. != "Missing")) %>% #remove any sites with missing data
      group_by_all() %>% 
      count() %>% #count of all patterns in rows
      ungroup()
    
    rm(withinNeutral)
    
    #Check whether there are any neutral sites at all, if there are not then print no useable sites
    if (dim(GTCounts)[1] > 1) {
      #Replace pattern with value
      GTCounts[GTCounts=="DerHom"]<- 2
      GTCounts[GTCounts=="Het"]<- 1
      GTCounts[GTCounts=="AncHom"]<- 0
      names(GTCounts)[7] = "Count"
      
      #count of the frequency of each pattern (which is the count of snps)
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
    } else{
      message1 = sprintf("No Useable Neutral Sites %s", chromo) #Chromosome 9 has nothing left
      print(message1)
    }
   
  }

  write.table(summaryInfo, file = paste(i, "_NeutralSFS_allChroms_SummaryFile_N6.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


}
