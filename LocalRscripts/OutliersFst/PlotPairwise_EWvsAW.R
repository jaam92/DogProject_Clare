library(tidyverse)
library(GenomicRanges)
library(boot)
library(ggrepel)

#read file in
pairwiseFst = read_delim("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/EW_vs_AW_allSites_rmNAN.weir.fst", delim = "\t") %>%
  mutate(POS_END = POS,
         updatedFst = ifelse(WEIR_AND_COCKERHAM_FST < 0, as.numeric(0), WEIR_AND_COCKERHAM_FST),
         bin = row_number()) 

genes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1.bed", sep = "\t", stringsAsFactors = F)
gene_names = read.table("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt", stringsAsFactors = F)

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

#Compute Fst within genes and figure out which window each position falls 
FstPosRanges = with(pairwiseFst, GRanges(CHROM, IRanges(start=POS, end = POS_END), score=updatedFst))

geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))

overlaps = findOverlaps(query = FstPosRanges, subject = geneRanges, type = "within") %>%
  as.data.frame()

#####Compute Fst per gene
SNPsPerGene = pairwiseFst %>%
  mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
  na.omit() %>% #remove variants that don't fall within genes
  group_by(range) %>%
  filter(n() > 20) %>%#remove genes with less than twenty snps
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
  ungroup() 

FstPerGene$SNPCount = SNPsPerGene$n[match(FstPerGene$range, SNPsPerGene$range)]

#Plot correlation between number of SNPs and Fst
CREBBP = FstPerGene %>% filter(GeneName =="CREBBP") %>% 
  select(FstNormSNPCount) %>% 
  as.numeric()#FST for CREBBP

FstPerGene$Label = ifelse(FstPerGene$FstNormSNPCount >= quantile(FstPerGene$FstNormSNPCount,prob=0.95), "Top5%","Remainder")

correlationFstSnpCount_CREBBP = ggplot(FstPerGene, aes(x=FstNormSNPCount,y=SNPCount,colour=Label)) + 
  geom_point() +  
  scale_colour_manual(name = "Status", values=c("blue","red")) + 
  geom_text_repel(size = 6, aes(label=ifelse((FstNormSNPCount == CREBBP) ,as.character(GeneName),'')), show.legend = F) + 
  theme_bw() +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size=24, face = "bold"), 
        plot.title = element_text(size=24, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=20))

#plot histogram of Fst and label EPAS 1
FstHistogram_CREBBP = ggplot() + 
  geom_histogram(data = FstPerGene, aes(x=FstNormSNPCount), bins=50) + 
  geom_vline(xintercept = CREBBP, colour="purple") + 
  theme_bw()

sum(FstPerGene$FstNormSNPCount >= CREBBP)/nrow(FstPerGene) #pvalue

####Compute Pi per gene EW and AW####
#Read file in 
AW_PI = read_delim("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/AW_allSites.sites.pi", delim = "\t") %>%
  mutate(bin = row_number()) 

piPosRanges = with(AW_PI, GRanges(CHROM, IRanges(start=POS, end = POS), score=PI))

geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))

overlaps = findOverlaps(query = piPosRanges, subject = geneRanges, type = "within") %>%
  as.data.frame()

#####Compute PI per gene for AW
AW_piPerGene = AW_PI %>%
  mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
  na.omit() %>% #remove variants that don't fall within genes
  group_by(range) %>%
  filter(n() > 20) %>%#remove genes with less than twenty snps
  summarise(sumPI = sum(as.numeric(PI)),
            PINormSNPCount = mean(as.numeric(PI))) %>% 
  mutate(PINormTranscriptLen = sumPI/(GeneSet$transcript_length[match(range,GeneSet$bin)]),
         GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) %>% #divide by transcript length
  ungroup()

#Read file in 
EW_PI = read_delim("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/EW_allSites.sites.pi", delim = "\t") %>%
  mutate(bin = row_number()) 

piPosRanges = with(EW_PI, GRanges(CHROM, IRanges(start=POS, end = POS), score=PI))

geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))

overlaps = findOverlaps(query = piPosRanges, subject = geneRanges, type = "within") %>%
  as.data.frame()

#####Compute PI per gene for EW
EW_piPerGene = EW_PI %>%
  mutate(range = overlaps$subjectHits[match(bin,overlaps$queryHits)]) %>% #column with corresponding gene from gene set
  na.omit() %>% #remove variants that don't fall within genes
  group_by(range) %>%
  filter(n() > 20) %>%#remove genes with less than twenty snps
  summarise(sumPI = sum(as.numeric(PI)),
            PINormSNPCount = mean(as.numeric(PI))) %>% 
  mutate(PINormTranscriptLen = sumPI/(GeneSet$transcript_length[match(range,GeneSet$bin)]),
         GeneName = GeneSet$AbbrevName[match(range,GeneSet$bin)]) %>% #divide by transcript length
  ungroup()

#delete large input files
rm(AW_PI)
rm(EW_PI)

#Plot Pi per gene
CREBBP_AW = AW_piPerGene %>% 
  filter(GeneName =="CREBBP") %>% 
  select(PINormSNPCount) %>% 
  as.numeric()

CREBBP_EW = EW_piPerGene %>% 
  filter(GeneName =="CREBBP") %>% 
  select(PINormSNPCount) %>% 
  as.numeric()

plotAW_piPerGene_CREBBP = ggplot() + 
  geom_histogram(data = AW_piPerGene, aes(x=PINormSNPCount), bins=50) + 
  geom_vline(xintercept = CREBBP_AW, colour="purple") + 
  theme_bw()

plotEW_piPerGene_CREBBP = ggplot() + 
  geom_histogram(data = EW_piPerGene, aes(x=PINormSNPCount), bins=50) +
  geom_vline(xintercept = CREBBP_EW, colour="purple") + 
  theme_bw()

sum(AW_piPerGene$PINormSNPCount <= CREBBP_AW)/nrow(AW_piPerGene)
sum(EW_piPerGene$PINormSNPCount <= CREBBP_EW)/nrow(EW_piPerGene)

#####BOOTSTRAPPING FST #####

#Look at whether data is normally distributed, it isn't.
#qqnorm(FstPerGene$FstNormSNPCount)
#qqline(FstPerGene$FstNormSNPCount)

