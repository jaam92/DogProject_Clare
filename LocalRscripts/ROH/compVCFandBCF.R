#Load libraries
library(data.table)
library(ggplot2)
library(hablar)
library(dplyr)
library(mgsub)
library(cowplot)
#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/ROH")
Individuals = read.table("~/Documents/DogProject_Clare/LocalRscripts/Dogs2Keep.txt")
nonQCROH_allCanids_bcfTools = read.delim("allPops_bcfTools_withGenMap.LROH")
#Load vcftools outputs
fnames = list.files(pattern = "\\ROHQC.LROH$")
df = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE, idcol = "FileName")#Read all the files and create a FileName column to store filenames
nonQCROH_allCanids_vcfTools = subset(df, df$CHROM != "CHROM") %>% retype() %>% as.data.frame()#remove header lines and convert character columns to numeric
rm(df)#delete old df

ggplot(nonQCROH_allCanids_bcfTools, aes(x=Length)) +  geom_histogram(bins = 200) + theme_bw() + xlim(0,1.50e6) + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18))

ggplot(nonQCROH_allCanids_vcfTools, aes(x=ROH_length)) +  geom_histogram(bins = 200) + theme_bw() + xlim(0,1.50e6)+ theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18))


#Split into classes
#Split into length intervals, find length for each indiv, add column in Mb, add range column
aggIntROH_vcfTools = nonQCROH_allCanids_vcfTools %>% select(ROH_length, INDV) %>% filter(ROH_length >= 100000 & ROH_length < 1000000) %>% group_by(INDV) %>% summarise(totalLen = sum(ROH_length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[0.1-1)Mb") %>% as.data.frame()

aggLongROH_vcfTools = nonQCROH_allCanids_vcfTools %>% select(ROH_length, INDV) %>% filter(ROH_length >= 1000000 & ROH_length < 10000000) %>% group_by(INDV) %>% summarise(totalLen = sum(ROH_length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[1-10)Mb")  %>% as.data.frame()

aggVLongROH_vcfTools = nonQCROH_allCanids_vcfTools %>% select(ROH_length, INDV) %>% filter(ROH_length >= 10000000) %>% group_by(INDV) %>% summarise(totalLen = sum(ROH_length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[10-63)Mb") %>% as.data.frame()

aggIntROH_bcfTools = nonQCROH_allCanids_bcfTools %>% select(Length, Sample) %>% filter(Length >= 100000 & Length < 1000000) %>% group_by(Sample) %>% summarise(totalLen = sum(Length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[0.1-1)Mb") %>% as.data.frame()

aggLongROH_bcfTools = nonQCROH_allCanids_bcfTools %>% select(Length, Sample) %>% filter(Length >= 1000000 & Length < 10000000) %>% group_by(Sample) %>% summarise(totalLen = sum(Length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[1-10)Mb")  %>% as.data.frame()

aggVLongROH_bcfTools = nonQCROH_allCanids_bcfTools %>% select(Length, Sample) %>% filter(Length >= 10000000) %>% group_by(Sample) %>% summarise(totalLen = sum(Length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[10-63)Mb") %>% as.data.frame()

#Merged interval data frames
FinalDF_vcfTools = rbind(aggIntROH_vcfTools, aggLongROH_vcfTools, aggVLongROH_vcfTools) %>% plyr::rename(c("INDV" = "Sample"))
FinalDF_bcfTools = rbind(aggIntROH_bcfTools, aggLongROH_bcfTools, aggVLongROH_bcfTools)

FinalDF_vcfTools$Range = factor(FinalDF_vcfTools$Range, levels=c("[10-63)Mb","[1-10)Mb","[0.1-1)Mb"))
FinalDF_vcfTools$Population = substr(FinalDF_vcfTools$Sample, 1,2)
FinalDF_vcfTools$Population = mgsub_dict(FinalDF_vcfTools$Population,conversions =  list("BC" = "Border Collie", "LB" = "Labrador Retriever", "PG" = "Pug", "TM" = "Tibetan Mastiff", "AW" = "Arctic Wolf", "EW" = "Ethiopian Wolf", "IR" = "Isle Royale"))
orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf", "Isle Royale")
FinalDF_vcfTools$Population = factor(FinalDF_vcfTools$Population, levels = orderPops)
FinalDF_vcfTools$Sample = factor(FinalDF_vcfTools$Sample, levels = Individuals$V1)

FinalDF_bcfTools$Range = factor(FinalDF_bcfTools$Range, levels=c("[10-63)Mb","[1-10)Mb","[0.1-1)Mb"))
FinalDF_bcfTools$Population = substr(FinalDF_bcfTools$Sample, 1,2)
FinalDF_bcfTools$Population = mgsub_dict(FinalDF_bcfTools$Population,conversions =  list("BC" = "Border Collie", "LB" = "Labrador Retriever", "PG" = "Pug", "TM" = "Tibetan Mastiff", "AW" = "Arctic Wolf", "EW" = "Ethiopian Wolf", "IR" = "Isle Royale"))
FinalDF_bcfTools$Population = factor(FinalDF_bcfTools$Population, levels = orderPops)
FinalDF_bcfTools$Sample = factor(FinalDF_bcfTools$Sample, levels = Individuals$V1)

#Plot as individual
vcf = ggplot() + geom_bar(data = FinalDF_vcfTools, aes(x=Sample, y=totalLen/10^9, fill=Range), stat="identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("[0.1-1)Mb"= "yellow", "[1-10)Mb" = "orange", "[10-63)Mb" = "red"), breaks = c("[0.1-1)Mb","[1-10)Mb", "[10-63)Mb"), name = "Range") + ylab("ROH Length per Bin (Gb)") + xlab("Group") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18))

bcf = ggplot() + geom_bar(data = FinalDF_bcfTools, aes(x=Sample, y=totalLen/10^9, fill=Range), stat="identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("[0.1-1)Mb"= "yellow", "[1-10)Mb" = "green", "[10-63)Mb" = "blue"), breaks = c("[0.1-1)Mb","[1-10)Mb", "[10-63)Mb"), name = "Range") + ylab("ROH Length per Bin (Gb)") + xlab("Group") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18))

plot_grid(vcf,bcf)

#Compare with ROH of at least 100,000 bb
compareMethods = merge(FinalDF_bcfTools, FinalDF_vcfTools, by = c("Sample", "Range"), all = TRUE)
compareMethods$Diff = compareMethods$totalLenGb.x - compareMethods$totalLenGb.y
table(ifelse(compareMethods$Diff > 0, "largerBCF", "largerVCF"))

IndividualLevel_vcfTools = compareMethods %>% 
  group_by(Sample) %>% 
  summarise(TotalROH = sum(as.numeric(totalLen.y), na.rm=TRUE)) %>%
  mutate(TotalROHMb = TotalROH/10^6, FROH = TotalROH/2500000000, Population = substr(Sample,1,2), Method = "vcfTools")

IndividualLevel_bcfTools = compareMethods %>%
  group_by(Sample) %>% 
  summarise(TotalROH = sum(as.numeric(totalLen.x), na.rm=TRUE)) %>%
  mutate(TotalROHMb = TotalROH/10^6, FROH = TotalROH/2500000000, Population = substr(Sample,1,2), Method = "bcfTools")

plot = rbind.data.frame(IndividualLevel_bcfTools, IndividualLevel_vcfTools)
plot$Population = mgsub_dict(plot$Population,conversions =  list("BC" = "Border Collie", "LB" = "Labrador Retriever", "PG" = "Pug", "TM" = "Tibetan Mastiff", "AW" = "Arctic Wolf", "EW" = "Ethiopian Wolf", "IR" = "Isle Royale"))
plot$Population = factor(plot$Population, levels = orderPops)

cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

dodge = position_dodge(width = 0.9)
ggplot(plot, aes(x=Population, y=FROH, colour = Method)) + theme_bw()  + scale_colour_manual(name = "Group", values = c(bcfTools="mediumpurple4", vcfTools="cyan4")) + geom_violin(size=1, position = dodge)  + theme_bw() + labs(x = "Group", y=expression(F[ROH])) + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18))




