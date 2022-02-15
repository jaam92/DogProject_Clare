#Load libraries
library(tidyverse)
library(GenomicRanges)

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

####Functions for Pi
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
    mutate(Population = popAbbrev) %>%
    ungroup()
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
    mutate(Population = popAbbrev) %>%
    ungroup()
  return(binData) 
}

####Functions for Fst
#Function to make dataframes and compute Fst
compFST = function(inFile){
  #read file in 
  pairwiseFst = read_delim(file = inFile, delim = "\t") %>%
    mutate(POS_END = POS,
           updatedFst = ifelse(WEIR_AND_COCKERHAM_FST < 0, as.numeric(0), WEIR_AND_COCKERHAM_FST),
           bin = row_number()) 
  
  #Compute Fst within genes and figure out which window each position falls 
  FstPosRanges = with(pairwiseFst, GRanges(CHROM, IRanges(start=POS, end = POS_END), score=updatedFst))
  geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))
  overlaps = findOverlaps(query = FstPosRanges, subject = geneRanges, type = "within") %>%
    as.data.frame()
  
  #Compute Fst per gene
  SNPsPerGene = pairwiseFst %>%
    mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
    na.omit() %>% #remove variants that don't fall within genes
    group_by(range) %>%
    filter(n() > 20) %>% #remove genes with less than twenty snps
    count() %>%
    mutate(GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) %>%
    ungroup()
  
  FstPerGene = pairwiseFst %>%
    mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
    na.omit() %>% #remove variants that don't fall within genes
    group_by(range) %>%
    filter(n() > 20) %>%#remove genes with less than twenty snps
    summarise(sumFst = sum(as.numeric(updatedFst)),
              FstNormSNPCount = mean(as.numeric(updatedFst))) %>% 
    mutate(FstNormTranscriptLen = sumFst/(GeneSet$transcript_length[match(range,GeneSet$bin)]),
           GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) %>% #divide by transcript length
    ungroup() %>%
    mutate(SNPCount = SNPsPerGene$n[match(range, SNPsPerGene$range)])
  
  return(FstPerGene)
}

#####Functions to make dataframes and compare fixed sites
compFixedSites = function(inFile){
  
  #read file in 
  fixedSites = read_delim(file = inFile, delim = "\t") %>%
    mutate(bin = row_number())
  
  #overlap genes
  fixedPosRanges = with(fixedSites, GRanges(CHROM, IRanges(start=POS, end = POS)))
  geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))
  overlaps = findOverlaps(query = fixedPosRanges, subject = geneRanges, type = "within") %>%
    as.data.frame()
  
  fixedSitesPerGene = fixedSites %>%
    mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
    na.omit() %>% #remove variants that don't fall within genes
    mutate(GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) 
  
  #count number of sites that are fixed
  CountPerGene = fixedSitesPerGene %>%
    group_by(GeneName) %>%
    count(name = "FixedSites") %>%
    ungroup()
  
  return(CountPerGene)
}

#####Functions for plotting summary stat around genes

#Fixed sites for 1Mb around any gene of interest
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
