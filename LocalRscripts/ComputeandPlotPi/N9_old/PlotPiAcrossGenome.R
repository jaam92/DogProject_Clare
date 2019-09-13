#Load Libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(data.table)
library(magrittr)
library(mgsub)

#Read in Files and a single data frame
setwd("~/Documents/DogProject_Clare/LocalRscripts/plotPi")
Fnames = lapply(Sys.glob("*downsampledN9_allChroms_genomeWide.pi"), read.table)
df = rbindlist(Fnames) %>% mutate(CompPi = as.numeric(V1)/as.numeric(V2), Population = substr(V3, 1, 2)) %>% set_colnames(c("totalPi", "totalSites", "Filename", "piPerChrom", "Population"))

#prep data for plotting
df$Population = mgsub_dict(df$Population,conversions =  list("BC" = "Border Collie", "LB" = "Labrador Retriever", "PG" = "Pug", "TM" = "Tibetan Mastiff", "AW" = "Arctic Wolf", "EW" = "Ethiopian Wolf", "IR" = "Isle Royale"))

orderPops = c("Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale")

df$Population = factor(df$Population, levels = orderPops)

cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")


#Regular Boxplot
dodge = position_dodge(width = 0.9)
(plotPi = ggplot(df, aes(x=Population, y= piPerChrom, colour = Population)) + geom_boxplot(size = 1.5, position = dodge) + geom_point(size=0.5) + scale_colour_manual(name = "Group", values = cbPalette) + theme_bw() + labs(x= "Group",y = expression(pi)) + ggtitle("Genome Wide Pi") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)))

#Average Pi across all chroms
df %>% group_by(Population) %>% summarise(mean = mean(piPerChrom))
