#Load libraries
library(data.table)
library(tidyverse)
library(mgsub)
library(RColorBrewer)

#Load in non-qc'd files
setwd("~/Documents/DogProject_Clare/LocalRscripts/ROH")
orderPops = c("border collie", "labrador retriever", "pug", "tibetan mastiff", "Arctic wolf", "Ethiopian wolf", "Isle Royale wolf")
Individuals = read_delim("~/Documents/DogProject_Clare/LocalRscripts/Dogs2Keep.txt", delim="\t", col_names = c("Sample"))

bcfTools_map = read_delim("allPops_bcfTools_withGenMap.LROH", delim="\t") %>% 
  filter(Sample%in%Individuals$Sample) %>%
  rename(ROH_length = Length) %>%
  mutate(method = "BCFTools_genetic_map") %>%
  select(Sample, ROH_length, method)

bcfTools_Nomap = read_delim("allPops_bcfTools.LROH", delim="\t") %>% 
  filter(Sample%in%Individuals$Sample) %>%
  rename(ROH_length = Length) %>%
  mutate(method = "BCFTools_no_map") %>%
  select(Sample, ROH_length, method)

fnames = list.files(pattern = "\\ROHQC.LROH$")
vcfTools = rbindlist(sapply(fnames, read_delim, delim="\t", simplify = FALSE), use.names = TRUE, idcol = "FileName") %>% 
  rename(Sample = INDV) %>%
  filter(CHROM != "CHROM" & Sample%in%Individuals$Sample)  %>%
  mutate(method = "VCFTools") %>%
  select(Sample, ROH_length, method)

allComparisons = rbind.data.frame(vcfTools, bcfTools_map, bcfTools_Nomap) %>%
  mutate(ROH_length = as.numeric(ROH_length),
         Population = substr(Sample, 1,2))

#Split into length intervals, find length for each indiv, add column in Mb, add range column
aggIntROH = allComparisons %>% 
  filter(ROH_length >= 100000 & ROH_length < 1000000) %>% 
  group_by(Sample, method) %>% 
  summarise(totalLen = sum(ROH_length)) %>% 
  mutate(totalLenGb = totalLen/10^9, Range = "[0.1-1)Mb") %>%
  ungroup()

aggLongROH = allComparisons %>% 
  filter(ROH_length >= 1000000 & ROH_length < 10000000) %>% 
  group_by(Sample, method) %>% 
  summarise(totalLen = sum(ROH_length)) %>% 
  mutate(totalLenGb = totalLen/10^9, Range = "[1-10)Mb") %>%
  ungroup()

aggVLongROH = allComparisons %>% 
  filter(ROH_length >= 10000000) %>% 
  group_by(Sample, method) %>% 
  summarise(totalLen = sum(ROH_length)) %>% 
  mutate(totalLenGb = totalLen/10^9, Range = "[10-63)Mb") %>%
  ungroup()

#Merged interval data frame
FinalDF = rbind.data.frame(aggIntROH, aggLongROH, aggVLongROH) %>% 
  mutate(Range = factor(Range, levels=c("[10-63)Mb","[1-10)Mb","[0.1-1)Mb")),
         Population = substr(Sample, 1,2),
         Population = str_replace_all(Population,pattern=c("BC" = "border collie", "LB" = "labrador retriever", "PG" = "pug", "TM" = "tibetan mastiff", "AW" = "Arctic wolf", "EW" = "Ethiopian wolf", "IR" = "Isle Royale wolf")),
         Population = factor(Population, levels = orderPops),
         method = gsub("_", " ", method)) 

##
compMethods = ggplot(FinalDF, aes(x=Sample, y=totalLenGb, fill=Range)) + 
  geom_col()  + 
  facet_wrap(~ method) +
  coord_flip() + 
  scale_fill_manual(values = c("[0.1-1)Mb"= "bisque3", 
                               "[1-10)Mb" = "darkgoldenrod",
                               "[10-63)Mb" = "indianred4"), 
                    breaks = c("[0.1-1)Mb","[1-10)Mb", "[10-63)Mb"), 
                    name = "Range") + 
  labs(y = "Total amount of genome in ROH(Gb)\nper length class", x = "Individual") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10), 
        plot.title = element_text(size=26, face = "bold", hjust=0.5),
        strip.text = element_text(size = 20), 
        axis.title = element_text(size=20),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=18))

#use color palette
cbPalette = c("Arctic wolf" = "gray25", "Ethiopian wolf" = "#D55E00",  "Isle Royale wolf" = "steelblue", "border collie" = "#009E73", "labrador retriever" = "gold3", "pug" = "mediumpurple4", "tibetan mastiff" = "#CC79A7")

FROH = rbind.data.frame(vcfTools, bcfTools_map, bcfTools_Nomap) %>%
  mutate(ROH_length = as.numeric(ROH_length)) %>%
  filter(ROH_length >= 10^6) %>%
  group_by(Sample, method) %>% 
  summarise_at(.vars = c("ROH_length"), sum, na.rm = TRUE) %>% 
  ungroup() %>%
  mutate(FROH = ROH_length/2500000000,
         Population = substr(Sample, 1,2),
         Population = str_replace_all(Population,pattern=c("BC" = "border collie", "LB" = "labrador retriever", "PG" = "pug", "TM" = "tibetan mastiff", "AW" = "Arctic wolf", "EW" = "Ethiopian wolf", "IR" = "Isle Royale wolf")),
         Population = factor(Population, levels = orderPops),
         method = gsub("_", " ", method))

dodge = position_dodge(width = 0.9)

compFROH = ggplot(FROH, aes(x=Population, y=FROH, colour = method)) + 
  scale_colour_brewer(name = "Method", palette="Set2") + 
  geom_violin(size=1, position = dodge)  + 
  geom_point(position = dodge, show.legend = F) +
  labs(x = "Population", y=expression(F[ROH])) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))


ggarrange(compMethods, compFROH, nrow = 2, ncol = 1)
