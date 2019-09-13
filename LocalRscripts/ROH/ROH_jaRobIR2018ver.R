#Load libraries
library(data.table)
library(ggplot2)
library(hablar)
library(dplyr)
library(mgsub)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/ROH")
Individuals = read.table("~/Documents/DogProject_Clare/LocalRscripts/Dogs2Keep.txt")

allCanidsTrueROH = read.delim("~/Documents/DogProject_Clare/LocalRscripts/ROH/jackieVersion_allChroms.LROH", col.names = c("CHROM", "AUTO_START","AUTO_END","MIN_START","MAX_END","N_VARIANTS_BETWEEN_MAX_BOUNDARIES","N_MISMATCHES","INDV")) %>%
  mutate(ROH_length = AUTO_END - AUTO_START)

short = allCanidsTrueROH %>% 
  filter(ROH_length >= 100000 & ROH_length < 1000000) %>% 
  group_by(INDV) %>%
  summarise_at(c("ROH_length"), sum, na.rm = T) %>%
  mutate(totalLenGb = ROH_length/10^9, 
         Range = "[0.1-1)Mb") %>% 
  as.data.frame()

medium = allCanidsTrueROH %>% 
  filter(ROH_length >= 1000000 & ROH_length < 10000000) %>% 
  group_by(INDV) %>%
  summarise_at(c("ROH_length"), sum, na.rm = T) %>% 
  mutate(totalLenGb = ROH_length/10^9, 
         Range = "[1-10)Mb")  %>% 
  as.data.frame()

long = allCanidsTrueROH %>% 
  filter(ROH_length >= 10000000) %>% 
  group_by(INDV) %>%
  summarise_at(c("ROH_length"), sum, na.rm = T) %>%
  mutate(totalLenGb = ROH_length/10^9, 
         Range = "[10-63)Mb") %>% 
  as.data.frame()

#Merge interval data frames
FinalDF = rbind(short, medium, long)
FinalDF$Range = factor(FinalDF$Range, levels=c("[10-63)Mb","[1-10)Mb","[0.1-1)Mb"))
FinalDF$Population = substr(FinalDF$INDV, 1,2)
FinalDF$Population = mgsub(FinalDF$Population, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))
FinalDF$Population = factor(FinalDF$Population, levels = orderPops)
FinalDF$INDV = factor(FinalDF$INDV, levels = Individuals$V1)

#Plot ROH ranges with all ROH greater than 100Kb
jackieROHRangePlot = ggplot(FinalDF, aes(x=INDV, y=totalLenGb, fill=Range)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  coord_flip() + 
  scale_fill_manual(values = c("[0.1-1)Mb"= "bisque3", 
                               "[1-10)Mb" = "darkgoldenrod",
                               "[10-63)Mb" = "indianred4"), 
                    breaks = c("[0.1-1)Mb",
                               "[1-10)Mb", 
                               "[10-63)Mb"), 
                    name = "Range") + 
  labs(x ="Individual", y = "ROH Length per Bin (Gb)") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10),
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) + 
  ggtitle("JaRob version")

