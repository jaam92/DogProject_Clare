#Set Libraries and working directory
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)
library(cowplot)

#Plot the priors
set.seed(1)
samps = seq(1:100000)
Ne1 = runif(samps, min = 10, max = 10000)
Ne2 = runif(samps, min = 11000, max = 70000)
Tbot = runif(samps, min = 10, max = 10000)
sampPriors = cbind.data.frame(samps,Ne1,Ne2,Tbot)
plotNe1 = ggplot() + geom_density(aes(x=Ne1), colour="blue") + ggtitle(expression(paste(Ne[1]))) + xlab("unif(10,10000)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotNe2 = ggplot() + geom_density(aes(x=Ne2), colour="red") + ggtitle(expression(paste(Ne[2]))) + xlab("unif(11000,70000)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotTbot = ggplot() + geom_density(aes(x=Tbot), colour="green") + ggtitle(expression(paste(T[bottleneck]))) + xlab("unif(10,10000)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

(plotPriorsTogether = plot_grid(plotNe1,plotNe2+ylab(""),plotTbot+ylab(""),nrow = 1, align = 'h'))


#Load ABC results
setwd("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/NeCurr10to10kNeAnc11kto70kTbot10to6k/")

fnames = list.files(pattern = "\\.txt$")

#Read all the files and create a FileName column to store filenames
results = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE, idcol = "FileName")

#Add some column names
#headerLine = c('SGETaskID',"SimNeOne","SimNeTwo","SimTime",'DistMetric1','DisMetricPaper')
#colnames(results) <- headerLine
names(results)[1]="SGETaskID"
names(results)[2]="SimNeOne"
names(results)[3]="SimNeTwo"
names(results)[4]="SimTime"
names(results)[5]="DistMetric1"
names(results)[6]="DisMetricPaper"

#Check whether any of the sims have the same Ne1,Ne2,T
dim(distinct(results)) == dim(results) #They dont

#Get the lowest 100 values on V5 which was the distance metric from Jackie's paper
###FIGURE OUT WHAT TO DO WITH TIES####
Top100 = results %>%
  top_n(-100, DisMetricPaper) 

#Plotting Time
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
modePosteriorNe1 = getmode(Top100$SimNeOne)
modePosteriorNe2 = getmode(Top100$SimNeTwo)
modePosteriorTbot = getmode(Top100$SimTime)

plotPosteriorNe1 = ggplot() + geom_density(data=sampPriors, aes(x=Ne1), colour="blue",size=1) + geom_density(data=Top100, aes(x=SimNeOne), colour="black") + geom_vline(data = Top100, xintercept = modePosteriorNe1, colour="purple") + ggtitle(expression(paste(Ne[1]))) + xlab(paste("Mode of Posterior = ", modePosteriorNe1)) + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotPosteriorNe2 = ggplot() + geom_density(data=sampPriors, aes(x=Ne2), colour="red",size=1) + geom_density(data=Top100, aes(x=SimNeTwo), colour="black") + geom_vline(data = Top100, xintercept = modePosteriorNe2, colour="purple") + ggtitle(expression(paste(Ne[2]))) + xlab(paste("Mode of Posterior = ", modePosteriorNe2)) + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotPosteriorTbot = ggplot() + geom_density(data=sampPriors, aes(x=Tbot), colour="green",size=1) + geom_density(data=Top100, aes(x=SimTime), colour="black") + geom_vline(data = Top100, xintercept = modePosteriorTbot, colour="purple") + ggtitle(expression(paste(T[bottleneck]))) + xlab(paste("Mode of Posterior = ", modePosteriorTbot)) + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

(plotPosteriorTogether = plot_grid(plotPosteriorNe1 + ggtitle(""),plotPosteriorNe2+ ggtitle("") +ylab(""),plotPosteriorTbot+ ggtitle("") +ylab(""),nrow = 1, align = 'h'))

(plot_grid(plotPriorsTogether,plotPosteriorTogether, labels = "AUTO", ncol = 1, align = 'v', axis = 'r'))

HetBins = Top100 %>%
  select(starts_with('V')) %>% 
  mutate_all(funs(sum(as.numeric(as.character(.)), na.rm = TRUE))) %>%  
  distinct() %>%
  gather() %>%
  mutate(proportion = value / sum(value),
         bin = as.numeric(row_number()-1)) 

(HistogramPlotPosteriorHets = ggplot() + geom_bar(data=HetBins, aes(x=bin, y=proportion), stat = "identity") + theme_bw() + labs(x="Count of Heterozygotes per 1Kb", y="Proportion of Total Hets") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom"))
