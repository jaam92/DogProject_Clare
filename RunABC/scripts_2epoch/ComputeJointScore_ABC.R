suppressPackageStartupMessages(library(tidyverse))
library(dplyr)
library(optparse)
library(data.table)

#Input parameters from command line if not found set to default
option_list <- list( 
  make_option(c("--LengthThreshold"), type="integer", default=200, 
              help="Minimum length threshold for simulation [default %default]",
              metavar="number"),
  make_option("--SGETaskID", type="integer", 
              help="SGE Task ID for job array"),
  make_option(c("--EmpInfile"), type="character", default="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt", 
              help="empirical dataset file name", metavar="character"),
  make_option(c("--SimInfile"), type="character", default=NULL, 
              help="simulated dataset file name", metavar="character"),
  make_option(c("--outFilePath"), type="character", default="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/jointDis",
              help="specify output file path", metavar="character"),
  make_option(c("--paramFile"), type="character", 
              help="simulation parameter file name must include SGE task ID", metavar="character")
)

#Parse the options
opt = parse_args(OptionParser(option_list=option_list))

#Input parameters
LengthThreshold = opt$LengthThreshold
SGETaskID = opt$SGETaskID
outFilePath = opt$outFilePath
empiricalData = fread(file=opt$EmpInfile) %>% 
  filter(sites_passing > LengthThreshold) 
simulatedData = read.table(file=opt$SimInfile, col.names = c("Pi","windowedS"))
simulationParameterFile = read.table(file=opt$paramFile, col.names = c("NeCurr","NeAnc","Tbot"))

#Function to compute score per bin
computeBinScore <- function(empiricalCol, simulatedCol){
  binScore = abs(empiricalCol - simulatedCol) #Diego's version
  #binScore = ifelse(empiricalCol != 0 , abs(empiricalCol - simulatedCol)/(empiricalCol), 0)
  return(as.numeric(binScore))
}

#Function to reformat seg site count data 
reformatDataFrame <- function(dataFrame, grp.var) {
  grp.var = enquo(grp.var)
  dataFrame %>%
    group_by(!!grp.var) %>%
    count() %>% #we want counts of occurances of SegSites per replicate
    ungroup()
}

#Function to reformat summary stat output 
createSumStatDF <- function(dataFrame, colInterest, summaryID){
  dataFrame %>%
    select(bin, colInterest) %>%
    group_by_at(vars(-c(colInterest))) %>%  # group by everything other than the value column. 
    mutate(row_id=1:n()) %>% 
    ungroup() %>%  # build group index
    spread(bin, colInterest) %>%
    select(-c(row_id)) %>%
    mutate(SumStat = summaryID) %>%
    cbind.data.frame(simulationParameterFile, TrackSumStats) %>%
    select(SumStat, NeCurr, NeAnc, Tbot, jointScore, SegSiteScore, PiScore, everything())
}

#Reformat data
simulatedSegSites = reformatDataFrame(simulatedData, windowedS)
empiricalSegSites = reformatDataFrame(empiricalData, windowedS)
combinedSegSites = merge(empiricalSegSites, simulatedSegSites, by="windowedS", all= T) %>% 
  dplyr::rename(Emp = n.x, Sim = n.y) #combine data and rename cols
combinedSegSites[is.na(combinedSegSites)] <- 0 #replace NA with 0
combinedSegSites$binScore = computeBinScore(combinedSegSites$Emp, combinedSegSites$Sim) #compute bin score

#Compute windowed pi
empiricalData$range = seq.int(nrow(empiricalData)) #new range numbers based on data used
simulatedData$bin = seq.int(nrow(simulatedData)) #add bin number for simulated data
simulatedData$windowedPi = simulatedData$Pi/as.numeric(empiricalData$sites_passing[match(simulatedData$bin,empiricalData$range)])#divide pi by the length of simulated data

#Bin Pi starting from 0 to 0.02
#include lowest to include 0
simulatedData$PIranges = as.numeric(as.character(cut(simulatedData$windowedPi, seq(from = 0, to = 0.025, by = 0.0005), include.lowest = T, labels = c(0:49))))
empiricalData$PIranges = as.numeric(as.character(cut(empiricalData$windowedPi, seq(from = 0, to = 0.025, by = 0.0005), include.lowest = T, labels = c(0:49))))

#Reformat data and compute scores
simulatedPI = reformatDataFrame(simulatedData, PIranges)
empiricalPI = reformatDataFrame(empiricalData, PIranges)

combinedPI = merge(empiricalPI, simulatedPI, by="PIranges", all= T) %>%   dplyr::rename(Emp = n.x, Sim = n.y) #combine data and rename cols
combinedPI[is.na(combinedPI)] <- 0 #replace NA with 0
combinedPI$binScore = computeBinScore(combinedPI$Emp, combinedPI$Sim) #compute bin score 

#Compute pi and S for simulated data
TrackSumStats = simulatedData %>% 
  summarise_at(c("windowedS","windowedPi"), sum,na.rm=TRUE) %>% 
  dplyr::rename(TotalS = windowedS, TotalPi = windowedPi) %>%
  as.data.frame()

#Merge both scores together and make sure all bins from 0 to max are present
mergedPIandSegSites = merge(combinedSegSites, combinedPI, by.x= "windowedS", by.y = "PIranges", all = T) %>%
  dplyr::rename(bin = windowedS) %>%
  complete(bin = seq(min(bin), max(bin))) %>% #fill in missing bins
  mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .))) #replace NAs with 0

mergedPIandSegSites$jointBinScore = mergedPIandSegSites$binScore.x + mergedPIandSegSites$binScore.y #compute joint score

#Create output file
simulationParameterFile$jointScore = sum(mergedPIandSegSites$jointBinScore) #Add the joint score to the simulation parameters file
simulationParameterFile$SegSiteScore = sum( mergedPIandSegSites$binScore.x) #Add the SegSite score to the simulation parameters file
simulationParameterFile$PiScore = sum(mergedPIandSegSites$binScore.y) #Add the Pi score to the simulation parameters file

hetOutput = createSumStatDF(mergedPIandSegSites, "Sim.x", "SegSitesPerBin")
piOutput = createSumStatDF(mergedPIandSegSites, "Sim.y", "piPerBin")

OutputFile = rbind.data.frame(hetOutput, piOutput)

write.table(OutputFile, file = paste(outFilePath,"ABCOutput_NoError_", SGETaskID, "_multSampJointScore.txt", sep=""), append = T, quote = F, sep = "\t", col.names = F, row.names = F)
