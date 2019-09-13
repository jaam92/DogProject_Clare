#Set Libraries and working directory
library(data.table)
library(tidyverse)
library(cowplot)

FixParamsResults = read.delim("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/results_3epoch/ResultsABC_simulateUnderBestModel/ABCOutput_NoError_50_multSampJointScore.txt", fill = T, sep="\t", header = F) 

#Add some column names
names(FixParamsResults)[1] = 'SummaryStat'
names(FixParamsResults)[2] = 'SimNeOne'
names(FixParamsResults)[3] = 'SimNeTwo'
names(FixParamsResults)[4] = 'SimTimeTwo'
names(FixParamsResults)[5] = 'SimTimeOne'
names(FixParamsResults)[6] = 'SimNeThree'
names(FixParamsResults)[7] = 'jointScore'
names(FixParamsResults)[8] = 'SegSiteScore'
names(FixParamsResults)[9] = 'PiScore'


HetBins = Top100 %>%
  filter(SummaryStat == "SegSitesPerBin") %>%
  select(starts_with('V')) %>% 
  mutate_all(funs(sum(as.numeric(as.character(.)), na.rm = TRUE))) %>%
  distinct() %>%
  gather() %>%
  mutate(proportion = value / sum(value),
         bin = as.numeric(1:n()-1),
         Data = "Simulated") 

(HistogramPlotPosteriorHets = ggplot() + geom_bar(data=HetBins, aes(x=bin, y=proportion), stat = "identity") + theme_bw() + labs(x="Count of Segregating Sites per 1Kb", y="Proportion of Total Segregating Sites") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom"))

#Compare to real data
realData = read.delim("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt") %>% filter(sites_total >= 200)

HetsRealData = realData %>%
  group_by(windowedS) %>%
  count() %>%
  ungroup() %>%
  mutate(propSites = n/sum(n),
         bin = as.numeric(windowedS),
         Data = "Empirical")

(HistogramPlotRDHets = ggplot() + geom_bar(data=HetsRealData, aes(x=windowedS, y=propSites), stat = "identity") + theme_bw() + labs(x="Count of Segregating Sites per 1Kb", y="Proportion of Total Segregating Sites") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom"))

#Plot Together
plot_grid(HistogramPlotRDHets + xlim(-1,10) + ggtitle("Empirical"), HistogramPlotPosteriorHets + xlim(-1,10) + ggtitle("Simulations") + ylab(""))

#Combine data
subsetReal = HetsRealData %>% select("bin","propSites", "Data")
subsetSim = HetBins %>% select("bin","proportion", "Data")
matchCols = c("bin","propSites","Data")
colnames(subsetReal) <- matchCols
colnames(subsetSim) <- matchCols
comboDF = rbind.data.frame(subsetReal,subsetSim)

ggplot(data=comboDF, aes(x=bin, y=propSites,fill=Data)) + geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + theme_bw() + labs(x="Count of Segregating sites per 1Kb", y="Proportion of Total Segregating sites") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom") + scale_fill_manual(values = c("Empirical"="blue","Simulated"="red"))

ggplot(data=comboDF, aes(x=bin, y=propSites,fill=Data)) + geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + theme_bw() + labs(x="Count of Segregating sites per 1Kb", y="Proportion of Total Segregating sites") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom") + scale_fill_manual(values = c("Empirical"="blue","Simulated"="red")) + xlim(0.5,8.5) + ylim(0,0.02)
