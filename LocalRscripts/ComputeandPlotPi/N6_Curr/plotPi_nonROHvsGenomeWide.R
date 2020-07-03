#Load Libraries
library(tidyverse)
library(data.table)
library(mgsub)
library(reshape2)

#Grab input files
setwd("~/Documents/DogProject_Clare/LocalRscripts/ComputeandPlotPi/N6_Curr/")
fnames = list.files(pattern = "\\_SummaryFile_N6.txt$")

#Generate data frame
##create columns with fileName, population, and compute pi 
df = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FileName") %>% 
  mutate(Population = gsub("_SummaryFile_N6.txt", "",FileName),
         nonROH = as.numeric(PI_nonROH_allSites)/as.numeric(CountSites_nonROH),
         ROH = as.numeric(PI_ROH_allSites)/as.numeric(CountSites_ROH),
         Exons = as.numeric(PI_Exon_allSites)/as.numeric(CountSites_Exons),
         genomeWide = as.numeric(PI_genomeWide_allSites)/as.numeric(CountSites_allSites)) %>%
  as.data.frame()
df[is.na(df)] = 0 #replace nan with 0

#Examine differences in diversity genome wide vs sites outside roh
Differences = function(dataFrame, comp1, comp2){
  dataFrame %>% 
  Population_by(Population) %>% 
  summarise(mean1 = mean(comp1), 
            mean2 = mean(comp2)) %>% 
  mutate(LargerComp1 = ifelse(mean1 > mean2, "T", "F"),
         Difference = mean1 - mean2)
  }
#Plot to see whether to run t-test or wilcoxon rank sum on differences
##not normal distribution 
#for(i in 1:7){print(ggplot(splitDF[[i]], aes(x=value,colour=variable)) + geom_density())}

#Paired wilcoxon rank sum to check whether difference between two categories is significant
PairedWilcoxonTest = function(dataFrame, comp1, comp2){
  dataFrame %>%
  select(Population, comp1, comp2) %>%
  melt(id=c("Population")) %>%
  split(.$Population) %>%
  map(~ wilcox.test(value ~ variable, paired=T, data = .)) %>%
  map_dbl("p.value")
}
#prep data for plotting
plotDF = df %>%
  select(Population, nonROH,ROH,Exons,genomeWide) %>%
  melt(id=c("Population"))

plotDF$Population = mgsub(plotDF$Population, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))

orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale")

plotDF$Population = factor(plotDF$Population, levels = orderPops)

cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

#Regular Boxplot
dodge = position_dodge(width = 0.9)
plotPi = ggplot(plotDF, aes(x=variable, y=value, colour = Population)) + geom_boxplot(size = 1, position = dodge) + scale_colour_manual(name = "Population", values = cbPalette) + theme_bw() + labs(x= "Population",y = expression(pi)) + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), plot.title=element_text(size=24, face = "bold", hjust=0.5), axis.title=element_text(size=24), legend.title=element_text(size=20), legend.text=element_text(size=20))

plotGW = ggplot(plotDF %>% 
                filter(variable == "genomeWide"), 
                aes(x=Population, y=value, colour = Population)) + 
  geom_point(size = 0, shape = 15) +
  geom_boxplot(size = 1, position = dodge, show.legend = FALSE) +
  scale_colour_manual(name = "Population", values = cbPalette) + 
  theme_bw() + 
  labs(x= "Population",y = expression(pi)) + 
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        plot.title=element_text(size=24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

plotGW_newLegend = plotGW + guides(color = guide_legend(override.aes = list(size = 10)))
