#Load Libraries
suppressPackageStartupMessages(library(GenomicRanges))
library(dplyr)
library(data.table)

#Make input file for ABC
df = fread("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/NeutralRegionsEW/10kWin_10kStep/EWSamps_allChroms_NeutralRegions_het_10000win_10000step_Header.txt")

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
#write.table(newDF, file = "/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile_sepIndividualsperLine.txt",sep = "\t", row.names = F, col.names = T, quote = F)

#####Make a data frame where we collapse individuals so each 1Kb window is merged across individuals into a single row 
#Divide by four to average across the 4 individuals
collapsedDF = newDF %>%
  group_by(chromo,window_start) %>%
  summarise_all(funs(sum)) %>%
  mutate(sites_total=sites_total/4, 
         sites_passing=sites_passing/4,
         calls=trunc(calls/4),
         range = row_number())

#####Compute pi for neutral regions 
pi = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/NeutralRegionsEW/allChroms.sites.pi") %>%
  mutate(bin = row_number()) #add column with row number

windowsDF = collapsedDF %>%
  select(chromo, window_start) %>%
  mutate(window_end = window_start + 10000, #add 10Kb since that was the step size of sliding windows
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
  summarise(sumPi = sum(as.numeric(PI))) %>% 
  mutate(windowedPi = sumPi/(collapsedDF$calls[match(range,collapsedDF$range)])) %>% #divide by callable sites in region
  ungroup()

#####Count Seg Sites for neutral regions 
segSites = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/NeutralRegionsEW/allChroms.snpden") %>%
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


####Add column with pi and S back to emipirical data file
df_SummaryStats = collapsedDF %>%
  mutate(windowedPi = piPerWindow$windowedPi[match(range,piPerWindow$range)],
         windowedS = segSitesPerWindow$windowedSegSites[match(range,segSitesPerWindow$range)]) 

df_SummaryStats[is.na(df_SummaryStats)] <- 0 #replace NA with 0


#####Write ABC empirical data input file
write.table(df_SummaryStats, file = "/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/EWSamps_allChroms_NeutralRegions_het_10000win_10000step_ABCinputFile.txt",sep = "\t", row.names = F, col.names = T, quote = F)

