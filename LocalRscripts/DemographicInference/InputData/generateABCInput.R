#Load Libraries
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(data.table)

#Make input file for ABC
df = fread("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_Header.txt")

#####Make a data frame where each individual's window is a separate row
newDF = data.frame()
for(i in c("EW4","EW5","EW6","EW7")){
  subsetDF = df %>% 
  select(chromo, window_start, sites_total, sites_passing, contains(i))
  colnames(subsetDF) <- c("chromo", "window_start", "sites_total", "sites_passing", "calls", "missing", "hets", "homRef", "homAlt")
  newDF = rbind.data.frame(newDF,subsetDF)
  print(dim(subsetDF))
  }
#Write a file where each individual's 1Kb windows are separated as rows
#write.table(newDF, file = "~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile_sepIndividualsperLine.txt",sep = "\t", row.names = F, col.names = T, quote = F)

#####Make a data frame where we collapse individuals so each 1Kb window is merged across individuals into a single row
collapsedDF = newDF %>%
  group_by(chromo,window_start) %>%
  summarise_all(funs(sum))


#####Compute pi for neutral regions 
pi = read.delim("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/InputData/allChroms.sites.pi") %>%
  mutate(bin = row_number()) #add column with row number

windowsDF = collapsedDF %>%
  select(chromo, window_start) %>%
  mutate(window_end = window_start + 1000, #add 1Kb since that was the step size of windows
         bin = row_number()) 

#Figure out which window each position falls 
piPosRanges = with(pi, GRanges(CHROM, IRanges(start=POS, end = POS), score=PI))
windowRanges = with(windowsDF, GRanges(chromo, IRanges(start=window_start, end = window_end)))
overlaps = findOverlaps(query = piPosRanges, subject = windowRanges, type = "within") %>%
  as.data.frame()

#####Compute pi per window
piPerWindow = pi %>%
  mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding 1kb window  from empirical data
  group_by(range) %>%
  summarise(windowedPi = sum(as.numeric(PI))/(1000)) %>% #divide by 1000 since that is win length
  ungroup()

#####Count Seg Sites for neutral regions 
segSites = read.delim("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/InputData/allChroms.snpden") %>%
  mutate(END_POS = POS + 10, #add 10 since that was window step size
         bin = row_number()) #add column with row number

#Figure out which window each seg site range falls in 
segSitesPosRanges = with(segSites, GRanges(CHROM, IRanges(start=POS, end = END_POS), score=SNPCOUNT))
overlaps = findOverlaps(query = segSitesPosRanges, subject = windowRanges, type = "within") %>%
  as.data.frame()

segSitesPerWindow = segSites %>%
  mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding 1kb window from empirical data
  group_by(range) %>%
  summarise(windowedSegSites = sum(as.numeric(SNPCOUNT))) %>% #divide by 1000 since that is win length
  ungroup() %>%
  na.omit() #remove seg site ranges that do not have a corresponding 1Kb window from empirical data


####Add column with pi and S back to window file
windowsDF_addPi = windowsDF %>%
  mutate(windowedPi = piPerWindow$windowedPi[match(bin,piPerWindow$range)],
         windowedPi = replace_na(windowedPi, 0)) 

#Add column to the collapsed DF
collapsedDF$bin = seq(nrow(collapsedDF))#add column with row number
collapsedDF$windowedPi = windowsDF_addPi$windowedPi[match(collapsedDF$bin,windowsDF_addPi$bin)]#add pi

write.table(collapsedDF, file = "~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt",sep = "\t", row.names = F, col.names = T, quote = F)
