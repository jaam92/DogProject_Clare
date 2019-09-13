library(tidyverse)
library(GenomicRanges)
library(boot)
library(ggrepel)

#read file in
pairwiseFst = read_delim("~/Documents/DogProject_Clare/LocalRscripts/Fst/EW_vs_TM_allSites_rmNAN.weir.fst", delim = "\t") %>%
  mutate(POS_END = POS,
         updatedFst = ifelse(WEIR_AND_COCKERHAM_FST < 0, as.numeric(0), WEIR_AND_COCKERHAM_FST),
         bin = row_number()) 

genes = read.delim("~/Documents/DogProject_Jaz/CaseControlROH/EnsemblGenes_CanFam3.1.bed", sep = "\t")
gene_names = read.table("~/Documents/DogProject_Jaz/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt")

####Make file with genes of interest
####Filter to only chromosomes 1-38
####keep only longest transcript (need to use distinct too bc some transcripts have equal length and we only want to keep one entry)

#gene set for Ensembl
GeneSet = genes %>%
  mutate(AbbrevName = gene_names$V2[match(name, gene_names$V1)], transcript_length = transcriptionEnd - transcriptionStart) %>%
  group_by(AbbrevName) %>% 
  filter(transcript_length == max(transcript_length) & as.numeric(chrom) <= 38) %>%
  distinct(AbbrevName,.keep_all= TRUE)  %>% 
  select(name, AbbrevName, chrom, transcriptionStart, transcriptionEnd, transcript_length) %>%
  as.data.frame()
GeneSet$bin = seq.int(nrow(GeneSet)) #add bin

#delete genes
rm(genes)

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

#Plot FST 
EPAS1 = FstPerGene %>% filter(GeneName =="EPAS1") %>% 
  select(FstNormSNPCount) %>% 
  as.numeric()#FST for EPAS1

#plot correlation between number of SNPs and Fst
FstPerGene$Label = ifelse(FstPerGene$FstNormSNPCount >= quantile(FstPerGene$FstNormSNPCount,prob=0.95), "Top5%","Remainder")

correlationFstSnpCount_EPAS1 = ggplot(FstPerGene, aes(x=FstNormSNPCount,y=SNPCount,colour=Label)) + 
  geom_point() + 
  scale_colour_manual(name = "Status", values=c("blue","red")) + 
  geom_text_repel(size = 6, aes(label=ifelse((FstNormSNPCount == EPAS1) ,as.character(GeneName),'')), show.legend = F) + 
  theme_bw() +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size=24, face = "bold"), 
        plot.title = element_text(size=24, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=20))

#plot histogram of Fst and label EPAS 1
FstHistogram_EPAS1 = ggplot() + 
  geom_histogram(data = FstPerGene, aes(x=FstNormSNPCount), bins=50) + 
  geom_vline(xintercept = EPAS1, colour="purple") + 
  theme_bw()

sum(FstPerGene$FstNormSNPCount >= EPAS1)/nrow(FstPerGene) #pvalue

####Compute Pi per gene EW and TM####
#Read file in 
TM_PI = read_delim("~/Documents/DogProject_Clare/LocalRscripts/Fst/TM_allSites.sites.pi", delim = "\t") %>%
  mutate(bin = row_number()) 

piPosRanges = with(TM_PI, GRanges(CHROM, IRanges(start=POS, end = POS), score=PI))

geneRanges = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))

overlaps = findOverlaps(query = piPosRanges, subject = geneRanges, type = "within") %>%
  as.data.frame()

#####Compute PI per gene for TM
TM_piPerGene = TM_PI %>%
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
EW_PI = read_delim("~/Documents/DogProject_Clare/LocalRscripts/Fst/EW_allSites.sites.pi", delim = "\t") %>%
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
rm(TM_PI)
rm(EW_PI)

#Plot PI in TM and EW
EPAS1_TM = TM_piPerGene %>% 
  filter(GeneName =="EPAS1") %>% 
  select(PINormSNPCount) %>% 
  as.numeric()

EPAS1_EW = EW_piPerGene %>% 
  filter(GeneName =="EPAS1") %>% 
  select(PINormSNPCount) %>% 
  as.numeric()

plotTM_piPerGene = ggplot() + 
  geom_histogram(data = TM_piPerGene, aes(x=PINormSNPCount), bins=50) + 
  geom_vline(xintercept = EPAS1_TM, colour="purple") + 
  theme_bw()

sum(TM_piPerGene$PINormSNPCount <= EPAS1_TM)/nrow(TM_piPerGene)

plotEW_piPerGene = ggplot() + 
  geom_histogram(data = EW_piPerGene, aes(x=PINormSNPCount), bins=50) + 
  geom_vline(xintercept = EPAS1_EW, colour="purple") + 
  theme_bw()

sum(EW_piPerGene$PINormSNPCount <= EPAS1_EW)/nrow(EW_piPerGene)

#####BOOTSTRAPPING FST #####

#Look at whether data is normally distributed, it isn't
#qqnorm(FstPerGene$FstNormSNPCount)
#qqline(FstPerGene$FstNormSNPCount)

#centerFstPerGene = FstPerGene$FstNormSNPCount - mean(FstPerGene$FstNormSNPCount) + EPAS1

#Make function to get the mean
#perGene_mean <- function(data,indices){
#  return(mean(data[indices]))
#}

#Bootstrap centered data
#set.seed(805)
#results = boot(centerFstPerGene, perGene_mean, 10000)
#hist(results$t)

#Compute p-value
#Step 1 calculate distance from observed
#ObservedMean = mean(FstPerGene$FstNormSNPCount)
#DistanceFromObserved = abs(ObservedMean - EPAS1)
#MoreExtremeEPAS1 = EPAS1 + DistanceFromObserved

#P-value is the probability of being than x units away from EPAS1
#(sum(results$t < ObservedMean) + sum(results$t > MoreExtremeEPAS1))/10000

#one sample wilcoxon
#wilcox.test(FstPerGene$FstNormSNPCount, mu = EPAS1)$p.value

