#Load libraries
#library(tidyverse)
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
countFixDerHomGene <- function(dataFrame, chromNum, winStart, winEnd, windowLength, geneName, pop, fixedAlleleCount){
  
  #read file in 
  fixedSites = dataFrame %>%
    filter(!!sym(pop) == fixedAlleleCount & CHROM == chromNum) %>%
    mutate(bin = row_number())
  
  #overlap genes
  fixedPosRanges = with(fixedSites, GRanges(CHROM, IRanges(start=POS, end = POS)))
  geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))
  overlaps = findOverlaps(query = fixedPosRanges, subject = geneRanges, type = "within") %>%
    as.data.frame()
  
  fixedDerHomSitesPerGeneEW = fixedSites %>%
    mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
    na.omit() %>% #remove variants that don't fall within genes
    mutate(GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) 
  
  CountPerGene = fixedDerHomSitesPerGeneEW %>%
    filter(POS >= winStart &  POS <= winEnd) %>%
    group_by(GeneName) %>%
    count() %>%
    ungroup() %>%
    mutate(Label = ifelse(GeneName == geneName, T, F))
  
  plotGene = ggplot(CountPerGene, aes(x=GeneName, y=n, colour=Label)) + 
    geom_point(size = 2.5) +
    scale_colour_manual(values=c("blue","red"), guide = "none") + 
    geom_text_repel(size = 16, aes(label=ifelse((GeneName == geneName) ,as.character(GeneName),'')), show.legend = F) + 
    theme_bw() + 
    #ggtitle(paste("Genes within", windowLength, "of", geneName, sep = " ")) +
    labs(x = "Gene name", 
         y="Count of fixed derived homozygotes\n(per gene)") +
    theme(axis.text.x = element_text(hjust = 0.5, angle = 90, vjust = 0.8, size = 40), 
          axis.text.y = element_text(size = 40), 
          #plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
          axis.title = element_text(size = 42),
          strip.text = element_text(size = 42),
          legend.position = "none")
  return(plotGene) 
}