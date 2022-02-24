#Set Libraries and working directory
library(data.table)
library(tidyverse)
library(ggpubr)
library(glue)

#Load ABC results
results = read.table("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/results_3epoch/ResultsABC_3epoch_NAnc200k_200Ksims/allReplicates_multSampJointScore_3epoch_fillCols.txt", sep ="")

#Add some column names
names(results)[1] = 'Fname'
names(results)[2] = 'SummaryStat'
names(results)[3] = 'SimNeOne'
names(results)[4] = 'SimNeTwo'
names(results)[5] = 'SimTimeTwo'
names(results)[6] = 'SimTimeOne'
names(results)[7] = 'SimNeThree'
names(results)[8] = 'jointScore'
names(results)[9] = 'SegSiteScore'
names(results)[10] = 'PiScore'


Top200 = results %>%
  filter(SummaryStat == "SegSitesPerBin") %>%
  top_n(-200, jointScore) 

results$Status = ifelse(results$jointScore %in% Top200$jointScore, "Accepted", "Rejected")

#Plot it to make sure we did in fact select the lowest scores
# ggplot(results %>% filter(SegSiteScore < 300), 
#        aes(x=SegSiteScore, y=PiScore, colour = Status)) + 
#   geom_point() +
#   scale_colour_manual(values=c(Accepted = "blue", Rejected = "gray60")) + 
#   theme_bw()

#See if scores are correlated
fitScores = lm(data = results %>% filter(Status == "Accepted"), formula = SegSiteScore ~ PiScore)
summary(fitScores)

#Plot linear regression
#ggplot(results %>% filter(Status == "Accepted"), 
#       aes(x=SegSiteScore,y=PiScore)) +
#  geom_point() + 
#  geom_smooth(method = "lm") + 
#  theme_bw()

#Split results data frame by summary statistic
splitResults = split.data.frame(results, results$SummaryStat)

#Makes some new column names associated with bin number
replaceColnames = results %>% 
  select(starts_with("V")) %>% #pull cols of interest
  colnames() %>%
  as.data.frame()
names(replaceColnames)[1] = "OG"
replaceColnames$binNum = 0:(nrow(replaceColnames)-1) #add bin number
replaceColnames$binName = paste0("bin",0:(nrow(replaceColnames)-1)) #add new column name

S = splitResults$SegSitesPerBin %>%
  filter(SegSiteScore < 300) %>% #zoom in on xaxis
  select(starts_with("V"), "Status") %>%
  pivot_longer(!Status, names_to = "variable") %>% 
  mutate(binNum = replaceColnames$binNum[match(variable,replaceColnames$OG)],
         variable = replaceColnames$binName[match(variable,replaceColnames$OG)]) 

densityS = ggplot(S %>% filter(binNum <= 13),
                  aes(value, fill=Status)) + 
  geom_density() + 
  facet_wrap(~variable) +
  scale_fill_manual(values=c(Accepted="blue", Rejected="gray60")) +
  theme_bw() +
  ggtitle(expression(~italic(S))) +
  #ggtitle(expression(Theta[W])) +
  theme(plot.title = element_text(size=24, hjust = 0.5, face = "bold"))

P = splitResults$piPerBin %>%
  filter(SegSiteScore < 300) %>% #zoom in on xaxis
  select(starts_with("V"), "Status") %>%
  pivot_longer(!Status, names_to = "variable") %>% 
  mutate(binNum = replaceColnames$binNum[match(variable,replaceColnames$OG)],
         variable = replaceColnames$binName[match(variable,replaceColnames$OG)]) 

densityP = ggplot(P %>% filter(binNum <= 13),
       aes(value, fill=Status)) + 
  geom_density() + 
  facet_wrap(~variable) +
  scale_fill_manual(values=c(Accepted="blue", Rejected="gray60")) +
  theme_bw() +
  ggtitle(expression(pi)) +
  theme(plot.title = element_text(size=24, hjust = 0.5, face = "bold"))

ggarrange(densityS + labs(x="Count per Bin", y=""), 
          densityP + labs(x="Count per Bin", y=""), 
          ncol = 2, 
          nrow = 1, 
          common.legend = T, 
          legend = "bottom")

#Compare to real data
realData = read.delim("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/InputData/EWSamps_allChroms_NeutralRegions_het_1000win_1000step_ABCinputFile.txt") %>% 
  filter(sites_total >= 200)

HetsRealData = realData %>%
  group_by(windowedS) %>%
  count() %>%
  ungroup() %>%
  mutate(propSites = n/sum(n),
         bin = as.numeric(windowedS),
         Data = "Empirical")

PiRealData = realData %>%
  mutate(bin = as.numeric(as.character(cut(windowedPi, seq(from = 0, to = 0.025, by = 0.0005), include.lowest = T, labels = c(0:49))))) %>%
  group_by(bin) %>%
  count() %>%
  ungroup() %>%
  mutate(propSites = n/sum(n),
         Data = "Empirical")

#Make data frames for plotting
#create data frame with pi ranges 
ranges = cut(realData$windowedPi, seq(from = 0, to = 0.025, by = 0.0005), include.lowest = T)
bin = as.numeric(0:49)
rangesDF = as.character(levels(ranges)) %>% 
  cbind(bin) %>% 
  as.data.frame()
names(rangesDF)[1] = "ranges"

#Plot the empirical data
empiricalPlottingDF = HetsRealData %>% #want to use data frame with the most bins
  select(Data, bin, n) %>%
  mutate(valuePi = PiRealData$n[match(bin,PiRealData$bin)],
         valuePi = replace_na(valuePi, 0),
         piRanges = rangesDF$ranges[match(bin, rangesDF$bin)]) %>%
  dplyr::rename(valueS = n, binNum = bin) 

#Add everything to the seg sites data frame
simulationsPlot = S %>%
  mutate(Data = "Simulated",
         valuePi = P$value,
         piRanges = rangesDF$ranges[match(binNum, rangesDF$bin)]) %>%
  dplyr::rename(valueS = value) %>%
  filter(binNum %in% empiricalPlottingDF$binNum) #subset down to bins in real data 

#Labels for S and Pi ranges facet wrap 
empiricalPlottingDF$wrapLabels = glue('{empiricalPlottingDF$binNum}:{empiricalPlottingDF$piRanges}')
simulationsPlot$wrapLabels = glue('{simulationsPlot$binNum}:{simulationsPlot$piRanges}')

#Use gsub to parse the labels
ggplot() + 
  geom_point(data = simulationsPlot %>% filter(binNum < 13),  aes(x=valueS, y=valuePi, colour=Status, shape = Data)) + 
  geom_point(data = empiricalPlottingDF %>% filter(binNum < 13),  aes(x=valueS, y=valuePi, shape = Data), colour="red") +
  #facet_wrap(~variable) +
  #facet_wrap(~wrapLabels, scales = "free", labeller = label_bquote(Theta[W]==.(gsub(":.*", "",wrapLabels))~"&"~ pi==.(gsub(".*:", "",wrapLabels)))) +
  facet_wrap(~wrapLabels, scales = "free", labeller = label_bquote(italic(S)==.(gsub(":.*", "",wrapLabels))~"&"~ pi==.(gsub(".*:", "",wrapLabels)))) +
  scale_colour_manual(values=c(Accepted="blue", Rejected="gray60")) +
  scale_shape_manual(name = "Data", values = c(Empirical=8, Simulated = 16)) + 
  guides(color = guide_legend(order = 0),
         shape = guide_legend(order = 1)) +
theme_bw() +
  #labs(x=expression("Count per bin" ~ Theta[W]), y=expression("Count per bin" ~ pi)) +
  labs(x=expression("Count per bin" ~ italic(S)), y=expression("Count per bin" ~ pi)) +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size=24, face = "bold"), 
        plot.title = element_text(size=24, hjust = 0.5),
        strip.text = element_text(size=20),
        legend.title=element_blank(), 
        legend.text=element_text(size=20))
