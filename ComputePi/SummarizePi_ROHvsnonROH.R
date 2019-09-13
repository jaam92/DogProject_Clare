#Load Libraries
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(data.table)

#Load files
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ComputePi")
indivs = read.table("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/AW_indivList_downsample_N6.txt")
ROH = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/FinalROHQCd_min1Mb_allIndivs_Dec2018_DogProjClare.LROH")

#Empty data frame to fill with summary data
summaryInfo = data.frame()

for (chrom in 1:38){
  
  #Read file with pi per site in 
  sitesPi = fread(file = paste("AW_downsampledN6_chr", chrom ,".sites.pi", sep=""))
  
  #Generate query data frames
  sitesPi$POS_END = sitesPi$POS
  
  #Pull ROH for specific indviduals
  chromo = paste("chr",chrom, sep = "")
  ROH_AW = ROH %>%
    filter(INDV %in% indivs$V1 & CHROM == chromo) %>%
    select(CHROM, AUTO_START, AUTO_END) %>%
    distinct(AUTO_START, AUTO_END, .keep_all = TRUE)
  
  #Make GRanges object to find overlaps
  posRanges = with(sitesPi, GRanges(CHROM, IRanges(start=POS, end = POS_END)))
  rohRanges = with(ROH_AW, GRanges(CHROM, IRanges(start=AUTO_START, end = AUTO_END)))
  
  #Find Positions Overlapping ROH
  RemovePos = findOverlaps(query = posRanges, subject = rohRanges, type = "within") %>% 
    as.data.frame() %>% 
    distinct(queryHits) 
  
  #Keep Positions outside of ROH
  FinalPosition = sitesPi[-c(RemovePos$queryHits),] 
  
  #Pi for ROH and nonROH as well as count of sites 
  PI_nonROH = FinalPosition %>% 
    summarise(PI_nonROH_allSites = sum(PI)) 
  
  CountSites_nonROH = dim(FinalPosition)[1] 
  
  PI_genomeWide = sitesPi %>% 
    summarise(PI_genomeWide_allSites = sum(PI))
  
  CountSites_allSites = dim(sitesPi)[1]
  
  chromSummary = cbind(chrom, PI_nonROH, CountSites_nonROH, PI_genomeWide, CountSites_allSites)
  
  #Dataframe with all chromosomes
  summaryInfo = bind_rows(summaryInfo, chromSummary)
  
}
pop = unique(gsub('[[:digit:]]+', '', indivs$V1))

write.table(summaryInfo, file = paste(pop, "_SummaryFile_N6.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

