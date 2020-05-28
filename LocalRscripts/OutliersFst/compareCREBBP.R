#Load libraries
library(tidyverse)
library(GenomicRanges)
library(ggrepel)
library(ggpubr)
setwd("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst")

#gene set for Ensembl
genes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1.bed", sep = "\t")
gene_names = read.table("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt")

####Make file with genes of interest
####Filter to only chromosomes 1-38
####keep only longest transcript (need to use distinct too bc some transcripts have equal length and we only want to keep one entry)
GeneSet = genes %>% 
  mutate(chrom = gsub("chr", "", chrom), 
         AbbrevName = gene_names$V2[match(name, gene_names$V1)],
         transcript_length = transcriptionEnd - transcriptionStart) %>% 
  group_by(AbbrevName) %>% 
  filter(transcript_length == max(transcript_length) & chrom %in% 1:38) %>% 
  ungroup() %>%
  distinct(AbbrevName,.keep_all= TRUE)  %>% 
  select(name, AbbrevName, chrom, transcriptionStart, transcriptionEnd, transcript_length) %>%
  mutate(chrom = paste0("chr",chrom)) %>%
  as.data.frame()
GeneSet$bin = seq.int(nrow(GeneSet)) #add bin

rm(genes)#delete genes

#Function to make dataframes
makeDataFrames = function(inFile){
  #read file in 
  inFileDF = read_delim(file = inFile, delim = "\t") %>%
    mutate(bin = row_number()) 
  #Find overlaps with Granges
  geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))
  piPosRanges = with(inFileDF, GRanges(CHROM, IRanges(start=POS, end = POS), score=PI))
  overlaps = findOverlaps(query = piPosRanges, subject = geneRanges, type = "within") %>%
    as.data.frame()
  #compute PI per gene 
  piPerGene = inFileDF %>%
    mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
    na.omit() %>% #remove variants that don't fall within genes
    group_by(range) %>%
    filter(n() > 20) %>% #remove genes with less than twenty snps
    ungroup() %>%
    mutate(GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) 
  }

#Function to compute pi per gene
computePI = function(dataFrame, popAbbrev){
  piPerGene = dataFrame %>%
    group_by(range, GeneName) %>%
    summarise(meanPI = mean(as.numeric(PI))) %>% 
    mutate(Population = popAbbrev)
  return(piPerGene)
  }

#Function to make pi per snp in CREBBP
computePI_CREBBP = function(dataFrame, popAbbrev){
  piPerGene = dataFrame %>%
    filter(GeneName == "CREBBP")
  print(ggplot(piPerGene, aes(x=POS, y=PI)) + 
          geom_point() + 
          theme_bw() + 
          labs(x = "Position" ,y=expression(pi)) +
          ggtitle(popAbbrev)) #plot pi per snp in CREBBP
  binData = piPerGene %>% 
    mutate(binned_r_value = cut(PI, seq(from = 0, to = 1, by = 0.1), include.lowest = T)) %>%
    group_by(binned_r_value) %>%
    tally() %>% 
    mutate(Population = popAbbrev)
  return(binData) 
  }

#Generate Dataframes
dfBC = makeDataFrames("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/BC_allSites.sites.pi") %>%
  mutate(Population = "BC")
dfEW = makeDataFrames("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/EW_allSites.sites.pi") %>%
  mutate(Population = "EW")
dfTM = makeDataFrames("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/TM_allSites.sites.pi") %>%
  mutate(Population = "TM")
dfAW = makeDataFrames("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/AW_allSites.sites.pi") %>%
  mutate(Population = "AW")

#Compute Pi across all genes
BC_piPerGene = computePI(dfBC, "BC")
EW_piPerGene = computePI(dfEW, "EW")
TM_piPerGene = computePI(dfTM, "TM")
AW_piPerGene = computePI(dfAW, "AW")

#Finescale evaluation of Pi for CREBBP
BC_CREBBP = computePI_CREBBP(dfBC, "BC")
EW_CREBBP = computePI_CREBBP(dfEW, "EW")
TM_CREBBP = computePI_CREBBP(dfTM, "TM")
AW_CREBBP = computePI_CREBBP(dfAW, "AW")

mergedDF = rbind.data.frame(BC_piPerGene,TM_piPerGene,EW_piPerGene, AW_piPerGene) %>%
  mutate(chrom = gsub("chr", "", GeneSet$chrom[match(GeneName, GeneSet$AbbrevName)]))

#data frame with CREBBP
CREBBP = subset.data.frame(mergedDF, mergedDF$GeneName == "CREBBP") 

#plot as histogram
ggplot() +
  geom_density(data = AW_piPerGene, aes(x=meanPI), colour ="gray") +
  geom_density(data = TM_piPerGene, aes(x=meanPI), colour ="blue") +
  geom_density(data = BC_piPerGene, aes(x=meanPI), colour ="green") +
  geom_density(data = EW_piPerGene, aes(x=meanPI), colour ="red") +
  theme_bw() +
  xlab(expression("Mean" ~ pi ~ "per gene"))

ggplot() + 
  geom_histogram(data = mergedDF, aes(x=meanPI), bins=100) +
  geom_vline(data = CREBBP, aes(xintercept = meanPI), colour="purple") + 
  facet_wrap(~Population, scale = "free") +
  theme_bw() +
  xlab(expression("Mean" ~ pi ~ "per gene"))

#plot as dots
ggplot(CREBBP, aes(x=Population,y=meanPI)) + 
  geom_point() +  
  theme_bw() +
  ylab(expression("Mean" ~ pi ~ "CREBBP"))

allSitesDF = rbind.data.frame(dfBC, dfTM, dfEW, dfAW) %>%
  filter(GeneName == "CREBBP")

ggplot(allSitesDF, aes(x=POS, y=PI, colour=Population)) + 
  geom_point() + 
  theme_bw() + 
  labs(x = "Position" ,y=expression(pi)) 



####Compare fixed Sites
fixedSites = read_delim("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsAW_N9.txt", delim = "\t")
fixedSites$bin = seq.int(nrow(fixedSites)) #add bin

fixedPosRanges = with(fixedSites, GRanges(CHROM, IRanges(start=POS, end = POS)))

geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))

overlaps = findOverlaps(query = fixedPosRanges, subject = geneRanges, type = "within") %>%
  as.data.frame()

fixedSitesPerGene = fixedSites %>%
  mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
  na.omit() %>% #remove variants that don't fall within genes
  mutate(GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) %>%
  filter(GeneName%in%EW_piPerGene$GeneName)

CountPerGene = fixedSitesPerGene %>%
  group_by(GeneName) %>%
  count() %>%
  ungroup()

#Together
ggplot(CountPerGene, aes(x=n)) + 
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 212) + #Count for CREBBP
  theme_bw() +
  ggtitle("p-value = 0.049") +
  labs(x="Number of Fixed Sites Per Gene (EW vs AW)")


#Fixed DerHom EW
fixedDerHomSitesPerGeneEW = fixedSites %>%
  filter(EW == "18") %>%
  mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
  na.omit() %>% #remove variants that don't fall within genes
  mutate(GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) %>%
  filter(GeneName%in%EW_piPerGene$GeneName)

CountPerGene = fixedDerHomSitesPerGeneEW %>%
  group_by(GeneName) %>%
  count() %>%
  ungroup()

#Distribution across whole genome
ggplot(CountPerGene, aes(x=n)) + 
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 98) + #Count for CREBBP
  theme_bw() + 
  ggtitle("p-value = 0.09") +
  labs(x="Number of Derived Homozygous Fixed Sites Per Gene (EW)")


#Function to plot 1Mb around any gene of interest
countFixDerHomGene <- function(chromNum, winStart, winEnd, windowLength, geneName){
  #Look at genes within 1 Mb around gene of interest
  fixedDerHomSitesPerGeneEW = fixedSites %>%
    filter(EW == "18" & CHROM == chromNum) %>%
    mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
    na.omit() %>% #remove variants that don't fall within genes
    mutate(GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) %>%
    filter(GeneName%in%EW_piPerGene$GeneName)
  
  CountPerGene = fixedDerHomSitesPerGeneEW %>%
    filter(POS >= winStart &  POS <= winEnd) %>%
    group_by(GeneName) %>%
    count() %>%
    ungroup() %>%
    mutate(Label = ifelse(GeneName == geneName, T, F))
  
  plotGene = ggplot(CountPerGene, aes(x=GeneName, y=n, colour=Label)) + 
          geom_point() +
          scale_colour_manual(values=c("blue","red"), guide = FALSE) + 
          geom_text_repel(size = 6, aes(label=ifelse((GeneName == geneName) ,as.character(GeneName),'')), show.legend = F) + 
          theme_bw() + 
          ggtitle(paste("Genes within", windowLength, "of", geneName, sep = " ")) +
          labs(x = "Gene Name", 
               y="Number of Derived Homozygous Fixed Sites Per Gene (EW)") +
    theme(axis.text.x = element_text(size  = 14), 
          axis.text.y = element_text(size = 20), 
          axis.title = element_text(size=24, face = "bold"), 
          plot.title = element_text(size=24, hjust = 0.5), 
          legend.position = "none")
  return(plotGene) 
}

CREBBP_1Mb = countFixDerHomGene("chr6", 36412372, 38531822, "1Mb", "CREBBP")
GRB2_1Mb = countFixDerHomGene("chr9", 4109244, 6181472, "1Mb", "GRB2")
MKL1_1Mb = countFixDerHomGene("chr10", 23599593, 25752076, "1Mb", "MKL1")
PYGB_1Mb = countFixDerHomGene("chr23", 426984, 2491240, "1Mb", "PYGB")
CDK8_1Mb = countFixDerHomGene("chr25", 12058947, 14178073, "1Mb", "CDK8")

HIFOutliers = ggarrange(CREBBP_1Mb + labs(x = "", y = "", title = ""), 
              MKL1_1Mb + labs(x = "", y = "", title = ""),
              PYGB_1Mb + labs(x = "", y = "", title = ""),
              CDK8_1Mb + labs(x = "", y = "", title = ""),
              GRB2_1Mb + labs(x = "", y = "", title = "") + theme(axis.text.x = element_blank()),
              align = 'v',
              ncol = 2, 
              nrow = 3)

annotatedHIFOutliers = annotate_figure(HIFOutliers,
                left = text_grob("Number of Derived Homozygous Fixed Sites Per Gene", 
                                 size=20, 
                                 face = "bold",
                                 rot = 90),
                bottom = text_grob("Gene Name", 
                                   size=20, 
                                   face = "bold"))
#Save to pdf
#pdf(file = "~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/CompareHIFOutliers.pdf",height = 12, width = 36)
#print(annotatedHIFOutliers)
#dev.off()


#Function to compute pi per gene
piWindowedAboutGene = function(dataFrame, chromNum, winStart, winEnd, windowLength, geneName){
  
  piPerGene = dataFrame %>%
    filter(CHROM == chromNum & POS >= winStart &  POS <= winEnd) %>%
    group_by(range, GeneName) %>%
    summarise(meanPI = mean(as.numeric(PI))) %>%
    mutate(Label = ifelse(GeneName == geneName, T, F))
  
  plotGenePi = ggplot(piPerGene, aes(x=GeneName, y=meanPI, colour=Label)) + 
    geom_point() +
    scale_colour_manual(values=c("blue","red"), guide = FALSE) + 
    geom_text_repel(size = 6, aes(label=ifelse((GeneName == geneName) ,as.character(GeneName),'')), show.legend = F) + 
    theme_bw() + 
    ggtitle(paste("Genes within", windowLength, "of", geneName, sep = " ")) +
    labs(x = "Gene Name", 
         y=expression("Mean" ~ pi ~ "per gene")) +
    theme(axis.text.x = element_text(size  = 14), 
          axis.text.y = element_text(size = 20), 
          axis.title = element_text(size=24, face = "bold"), 
          plot.title = element_text(size=24, hjust = 0.5), 
          legend.position = "none")
  
  return(plotGenePi) 
}

CREBBP_1Mb_pi = piWindowedAboutGene(dfEW,"chr6", 36412372, 38531822, "1Mb", "CREBBP")
GRB2_1Mb_pi = piWindowedAboutGene(dfEW,"chr9", 4109244, 6181472, "1Mb", "GRB2")
MKL1_1Mb_pi = piWindowedAboutGene(dfEW,"chr10", 23599593, 25752076, "1Mb", "MKL1")
PYGB_1Mb_pi = piWindowedAboutGene(dfEW,"chr23", 426984, 2491240, "1Mb", "PYGB")
CDK8_1Mb_pi = piWindowedAboutGene(dfEW,"chr25", 12058947, 14178073, "1Mb", "CDK8")

