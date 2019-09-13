#Load libraries
library(data.table)
library(ggplot2)
library(hablar)
library(dplyr)
library(mgsub)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/ROH")
Individuals = read.table("~/Documents/DogProject_Clare/LocalRscripts/Dogs2Keep.txt")
df = read.delim("allPops_bcfTools.LROH") %>% retype() %>% mutate(Population = substr(Sample,1,2))

ggplot() + geom_density(data = df, aes(x=V8)) + labs(x="Quality Score") + ggtitle("Pre-Filtering") + theme_bw() + theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), axis.text.y = element_text(size =20), plot.title=element_text(size =24, face = "bold", hjust=0.5), axis.title=element_text(size=24))

#remove anything less than 1kb 
#aggregate by individual and population and find FROH

IndividualLevel = df %>% 
  filter(Length >= 1000) %>%
  group_by(Sample) %>% 
  summarise(TotalROH = sum(as.numeric(Length))) %>%
  mutate(TotalROHMb = TotalROH/10^6, FROH = TotalROH/2500000000, Population = substr(Sample,1,2))

IndividualsSampled = df %>%
  distinct(Sample) %>% 
  mutate(Population = substr(Sample, 1, 2)) %>%
  group_by(Population) %>%
  tally()

PopulationLevel = df %>%
  filter(Length >= 1000) %>%
  group_by(Population) %>% 
  summarise(TotalROH = sum(as.numeric(Length))) %>%
  mutate(TotalROHMb = TotalROH/10^6, FROH = TotalROH/(2500000000*IndividualsSampled$n[match(Population, IndividualsSampled$Population)]))

#Plot FROH
orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf", "Isle Royale")
IndividualLevel$Population = mgsub_dict(IndividualLevel$Population,conversions =  list("BC" = "Border Collie", "LB" = "Labrador Retriever", "PG" = "Pug", "TM" = "Tibetan Mastiff", "AW" = "Arctic Wolf", "EW" = "Ethiopian Wolf", "IR" = "Isle Royale"))
IndividualLevel$Population = factor(IndividualLevel$Population, levels = orderPops)

cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

ggplot(IndividualLevel, aes(x=Population, y=FROH, colour=Population)) + scale_colour_manual(name = "Group", values = cbPalette) + geom_boxplot(size=1) + geom_point(size=0.5) + theme_bw() + labs(x = "Group", y=expression(F[ROH])) + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18))

#Split into length intervals, find length for each indiv, add column in Gb, add range column
aggIntROH = df %>% select(Length, Sample) %>% filter(Length >= 100000 & Length < 1000000) %>% group_by(Sample) %>% summarise(totalLen = sum(Length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[0.1-1)Mb") %>% as.data.frame()

aggLongROH = df %>% select(Length, Sample) %>% filter(Length >= 1000000 & Length < 10000000) %>% group_by(Sample) %>% summarise(totalLen = sum(Length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[1-10)Mb")  %>% as.data.frame()

aggVLongROH = df %>% select(Length, Sample) %>% filter(Length >= 10000000) %>% group_by(Sample) %>% summarise(totalLen = sum(Length)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[10-63)Mb") %>% as.data.frame()

#Merged interval data frames
FinalDF = rbind(aggIntROH, aggLongROH, aggVLongROH)
FinalDF$Range = factor(FinalDF$Range, levels=c("[10-63)Mb","[1-10)Mb","[0.1-1)Mb"))
FinalDF$Population = substr(FinalDF$Sample, 1,2)
FinalDF$Population = mgsub_dict(FinalDF$Population,conversions =  list("BC" = "Border Collie", "LB" = "Labrador Retriever", "PG" = "Pug", "TM" = "Tibetan Mastiff", "AW" = "Arctic Wolf", "EW" = "Ethiopian Wolf", "IR" = "Isle Royale"))
FinalDF$Population = factor(FinalDF$Population, levels = orderPops)
FinalDF$Sample = factor(FinalDF$Sample, levels = Individuals$V1)

#Plot
ggplot(FinalDF, aes(x=Sample, y=totalLen/10^9, fill=Range)) + geom_bar(stat="identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("[0.1-1)Mb"= "yellow", "[10-63)Mb" = "blue", "[1-10)Mb" = "green"), breaks = c("[0.1-1)Mb","[1-10)Mb", "[10-63)Mb"), name = "Range") + ylab("ROH Length per Bin (Gb)") + xlab("Group") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18))
