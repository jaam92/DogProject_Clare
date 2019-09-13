#Load Libraries
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(data.table)

#set directory
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ComputePi")
#Load ROH
ROH = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/FinalROHQCd_min1Mb_allIndivs_Dec2018_DogProjClare.LROH")
#Load Exons and make 1 based
Exons = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/CanFam3.1_ExonIntronNeutral/Exons/EnsemblGenes_CanFam3.1_SingleTranscript.bed") %>%
  select(CHROM, exonStarts, exonEnds) %>%
  mutate(exonStarts1based = exonStarts+1)
#Load Neutral regions and make 1 based
Neutral = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/CanFam3.1_ExonIntronNeutral/NeutralRegions/allChroms_outside.ensemble_rm.phastConsEl_cutoff_0.4.bed") %>%
  mutate(START1based = START+1)
#Populations to iterate through
pops = c("BC","LB","PG","TM","AW","EW","IR")

#Fxn to compute pi and total sites in region
computePi = function(dataFrame, outputColnamePI, outputColnameCount){
  Pi = dataFrame %>%
    summarise(z = sum(PI))
  Count = dim(dataFrame)[1]
  combo = cbind(Pi,Count)
  colnames(combo) = c(outputColnamePI, outputColnameCount)
  return(combo)
}

for (i in pops){
  
  indivs = read.table(file = paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/" ,i, "_indivList_downsample_N6.txt", sep = ""))
  
  #Empty data frame to fill with summary data
  summaryInfo = data.frame()
  
  for (chrom in 1:38){
    
    #Read file with pi per site in 
    sitesPi = fread(file = paste(i,"_downsampledN6_chr", chrom ,".sites.pi", sep=""))
    
    #Generate query data frames
    sitesPi$POS_END = sitesPi$POS
    
    #Pull ROH for specific chrom and indviduals
    chromo = paste("chr",chrom, sep = "")
    subsetROH = ROH %>%
      filter(INDV %in% indivs$V1 & CHROM == chromo) %>%
      select(CHROM, AUTO_START, AUTO_END) %>%
      distinct(AUTO_START, AUTO_END, .keep_all = TRUE)
    #Pull exons for specific chrom
    subsetExons = Exons %>%
      filter(CHROM == chromo) %>%
      select(CHROM, exonStarts1based, exonEnds)
   #Pull neutral regions for specific chrom
    subsetNeutral = Neutral %>%
      filter(CHROM == chromo) %>%
      select(CHROM, START1based, END)
   
    #Make GRanges object to find overlaps
    posRanges = with(sitesPi, GRanges(CHROM, IRanges(start=POS, end = POS_END)))
    rohRanges = with(subsetROH, GRanges(CHROM, IRanges(start=AUTO_START, end = AUTO_END)))
    exonRanges = with(subsetExons, GRanges(CHROM, IRanges(start=exonStarts1based, end = exonEnds)))
    neutralRanges = with(subsetNeutral, GRanges(CHROM, IRanges(start=START1based, end = END)))

    #Find Positions Overlapping ROH
    POSROHOverlaps = findOverlaps(query = posRanges, subject = rohRanges, type = "within") %>% 
      as.data.frame() %>% 
      distinct(queryHits) 
    
    #Find Positions Overlapping Exons
    POSExonOverlaps = findOverlaps(query = posRanges, subject = exonRanges, type = "within") %>% 
      as.data.frame() %>% 
      distinct(queryHits)

    #Find Positions Overlapping Exons
    POSNeutralOverlaps = findOverlaps(query = posRanges, subject = neutralRanges, type = "within") %>%
      as.data.frame() %>%
      distinct(queryHits)
    
    #Keep Positions outside of ROH
    outsideROH = sitesPi[-c(POSROHOverlaps$queryHits),]
    
    #Keep Positions within ROH
    withinROH = sitesPi[c(POSROHOverlaps$queryHits),] 
      
    #Keep Positions within Exons
    withinExons = sitesPi[c(POSExonOverlaps$queryHits),]

    #Keep Positions within Neutral
    withinNeutral = sitesPi[c(POSNeutralOverlaps$queryHits),]
    
    #Pi for ROH and nonROH as well as count of sites 
    PIoutsideROH = computePi(outsideROH, "PI_nonROH_allSites", "CountSites_nonROH")
    PIwithinROH = computePi(withinROH, "PI_ROH_allSites", "CountSites_ROH")
    PIwithinExons = computePi(withinExons, "PI_Exon_allSites", "CountSites_Exons")
    PIwithinNeutral = computePi(withinNeutral, "PI_Neutral_allSites", "CountSites_Neutral")
    PIgenomeWide = computePi(sitesPi, "PI_genomeWide_allSites", "CountSites_allSites")
    
    chromSummary = cbind(chrom,PIoutsideROH,PIwithinROH,PIwithinExons,PIwithinNeutral,PIgenomeWide)
    
    #Dataframe with all chromosomes
    summaryInfo = bind_rows(summaryInfo, chromSummary)
    
  }

  write.table(summaryInfo, file = paste(i, "_SummaryFile_N6.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
}

