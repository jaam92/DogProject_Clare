#Set Libraries and working directory
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)

#Load ABC results
setwd("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/ABCResults_NeCurr900to10kNeAnc11kto40kT1700to2400/")
fnames = list.files(pattern = "\\.txt$")

#Read all the files and create a FileName column to store filenames
results = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE, idcol = "FileName")

#Add column names
headerLine = c('SGETaskID',"SimNeOne","SimNeTwo","SimTime",'DistMetric1','DisMetricPaper',"BIN0","BIN1","BIN2","BIN3","BIN4",'BIN5','BIN6','BIN7','BIN8','BIN9','BIN10','BIN11','BIN12','BIN13','BIN14','BIN15','BIN16','BIN17','BIN18','BIN19','BIN20','BIN21','BIN22','BIN23','BIN24','BIN25')
colnames(results) <- headerLine

#Check whether any of the sims have the same Ne1,Ne2,T
dim(distinct(results)) == dim(results) #They dont

#Get the lowest 100 values on V5 which was the distance metric from Jackie's paper
###FIGURE OUT WHAT TO DO WITH TIES####
Top100 = results %>%
  top_n(-100, DisMetricPaper) 

#Plotting Time
source(file = "~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/PriorsFromPSMC.R")
sampPriors = cbind.data.frame(samps,Ne1,Ne2,Tbot)

plotPosteriorNe1 = ggplot() + geom_density(data=sampPriors, aes(x=Ne1), colour="blue",size=1) + geom_density(data=Top100, aes(x=SimNeOne), colour="black") + ggtitle(expression(paste(Ne[1]))) + xlab("unif(900,10000)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotPosteriorNe2 = ggplot() + geom_density(data=sampPriors, aes(x=Ne2), colour="red",size=1) + geom_density(data=Top100, aes(x=SimNeTwo), colour="black") + ggtitle(expression(paste(Ne[2]))) + xlab("unif(11000,40000)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotPosteriorTbot = ggplot() + geom_density(data=sampPriors, aes(x=Tbot), colour="green",size=1) + geom_density(data=Top100, aes(x=SimTime), colour="black") + ggtitle(expression(paste(T[bottleneck]))) + xlab("unif(1700,2400)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotPosteriorTogether = plot_grid(plotPosteriorNe1 + ggtitle(""),plotPosteriorNe2+ ggtitle("") +ylab(""),plotPosteriorTbot+ ggtitle("") +ylab(""),nrow = 1, align = 'h')

plot_grid(plotPriorsTogether,plotPosteriorTogether, labels = "AUTO", ncol = 1, align = 'v', axis = 'r')

HetBins = Top100 %>%
  select(starts_with('BIN')) %>% 
  mutate_all(funs(sum(as.numeric(as.character(.)), na.rm = TRUE))) %>%  
  distinct() %>%
  gather() %>%
  mutate(proportion = value / sum(value),
         bin = as.numeric(as.character(gsub("BIN","", key)))) 

(HistogramPlotPosteriorHets = ggplot() + geom_bar(data=HetBins, aes(x=bin, y=proportion), stat = "identity") + theme_bw() + labs(x="Count of Heterozygotes per 1Kb", y="Proportion of Total Hets") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom"))
