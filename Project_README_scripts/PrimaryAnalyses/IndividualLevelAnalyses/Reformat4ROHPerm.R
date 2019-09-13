####This script will take in annotation data across all chromosomes for each individual then reformat the data so it can be used in Christian's script to permute ROH

####Input file:annotation file for each individual with data from all chroms
####Output files(3): 
	###annotation file for each individual with data from all chroms with snps as headers and their GT
	###annotation file for each individual with data from all chroms with snps as headers and whether site is within 10kb or greater ROH
	###annotation file for each individual with data from all chroms with snps as headers and whether site is within 1Mb or greater ROH

#Load Libraries
library(dplyr)
library(tidyr)

#Get input files
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT/AllChroms")
filenames = list.files(pattern = glob2rx("allIndivs_Chr*_annotatedGTwithVEP.txt"))

#Loop thru all individuals files
for (i in 1:length(filenames)){
  
  #Read files in 
  annots = read.delim(file = filenames[i])
  
  #Data to save
  indivGTInfo = annots %>%
    select(ID,POS,GT) %>%
    spread(POS,GT) %>%
    mutate(ID = gsub(".*annotatedGTwithVEP_\\s*|_Chr.*", "", ID, perl = T)) %>%
    arrange(ID)
  
  indivROHInfo = annots %>%
    select(ID,POS,withinROH_10Kb) %>%
    mutate(withinROH_10Kb = ifelse(withinROH_10Kb == 1, "ROH", "nonROH")) %>%
    spread(POS,withinROH_10Kb) %>%
    mutate(ID = gsub(".*annotatedGTwithVEP_\\s*|_Chr.*", "", ID, perl = T)) %>%
    arrange(ID)
  
  indivROHInfo_1Mb = annots %>%
    select(ID,POS,withinROH_1Mb) %>%
    mutate(withinROH_1Mb = ifelse(withinROH_1Mb == 1, "ROH", "nonROH")) %>%
    spread(POS,withinROH_1Mb) %>%
    mutate(ID = gsub(".*annotatedGTwithVEP_\\s*|_Chr.*", "", ID, perl = T)) %>%
    arrange(ID)

  indivROHInfo_smallPop = annots %>%
    select(ID,POS,withinROH_10Kb,withinROH_1Mb) %>% 
    mutate(NewIndicator = ifelse(withinROH_10Kb == "1" & withinROH_1Mb != "1", "ROH", "nonROH"),
           withinROH_10Kb = ifelse(withinROH_10Kb == 1, "ROH", "nonROH"),
           smallPopSize = ifelse(NewIndicator == "ROH" & withinROH_10Kb == "ROH", "ROH", "nonROH")) %>%
    select(-c(NewIndicator,withinROH_10Kb,withinROH_1Mb)) %>%   
    spread(POS,smallPopSize) %>%
    mutate(ID = gsub(".*annotatedGTwithVEP_\\s*|_Chr.*", "", ID, perl = T)) %>%
    arrange(ID)

  #Genotype data frame with all individuals and only segregating sites
  FilterIndivGT = Filter(function(x)(length(unique(x))>1), indivGTInfo)
  SegSites = colnames(FilterIndivGT)

  FilterIndivROH = indivROHInfo %>%
    select(SegSites)

  FilterIndivROH_1Mb = indivROHInfo_1Mb %>%
    select(SegSites)

 FilterIndivROH_smallPop = indivROHInfo_smallPop %>%
    select(SegSites)

  #Write to output files the only numbers in the file names should be 'ChrX' 
  #write.table(FilterIndivGT,  file= paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ROHPermutation/PermFormatGT_", filenames[i], sep=""), row.names = FALSE , col.names= TRUE, quote = FALSE, sep = "\t" )
  #write.table(FilterIndivROH,  file= paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ROHPermutation/PermFormat_TenKbROH_", filenames[i], sep=""), row.names = FALSE , col.names= TRUE, quote = FALSE, sep = "\t" )
  #write.table(FilterIndivROH_1Mb,  file= paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ROHPermutation/PermFormat_OneMbROH_", filenames[i], sep=""), row.names = FALSE , col.names= TRUE, quote = FALSE, sep = "\t" )
 write.table(FilterIndivROH_smallPop,  file= paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ROHPermutation/PermFormat_btwnTenKbandOneMbROH_", filenames[i], sep=""), row.names = FALSE , col.names= TRUE, quote = FALSE, sep = "\t" )
 
}

