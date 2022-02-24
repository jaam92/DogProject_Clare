#Set Libraries and working directory
library(tidyverse)
library(data.table)
library(hdrcde)
library(ggpubr)
options(scipen = 999)
#Plot the priors
set.seed(2)
samps = seq(1:100000)
Ne1 = runif(samps, min = 10, max = 1000)
Ne2 = runif(samps, min = 1000, max = 60000)
Ne3 = runif(samps, min = 40000, max = 200000)
Tbot = runif(samps, min = 10, max = 1000)
Tbot2 = runif(samps, min = 1500, max = 15000)
sampPriors = cbind.data.frame(samps,Ne1,Ne2,Ne3,Tbot,Tbot2)

#Plotting Function for priors
plotPrior <- function(xAxisTitle, title, parameter, color){
  ggplot() + 
    geom_density(aes(x=parameter), colour=color, outline.type = "full") + 
    ggtitle(title) + 
    xlab(xAxisTitle) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size  = 20), 
          axis.text.y = element_text(size = 20), 
          axis.title=element_text(size=20), 
          plot.title = element_text(size=20, hjust = 0.5))
}

plotNe1 = plotPrior("unif(10,1000)", expression(paste(Ne[CURRENT])), Ne1, "gray30")

plotNe2 = plotPrior("unif(1000,60000)", expression(paste(Ne[INTERMIDIATE])), Ne2, "gray30")

plotNe3 = plotPrior("unif(40000,200000)", expression(paste(Ne[ANCIENT])), Ne3, "gray30")

plotTbot = plotPrior("unif(10,1000)", expression(paste(TBOT[RECENT])), Tbot, "gray30")

plotTbot2 = plotPrior("unif(1500,15000)", expression(paste(TBOT[ANCIENT])), Tbot2, "gray30")

plotPriorsTogether = ggarrange(plotNe1, 
                               plotTbot + ylab(""), 
                               plotNe2 + ylab(""), 
                               plotTbot2 + ylab(""), 
                               plotNe3 + ylab(""), 
                               nrow = 1, 
                               ncol = 5, 
                               align = 'h')

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

#Grab the lowest 200 values 
Top200 = results %>%
  filter(SummaryStat == "SegSitesPerBin") %>%
  top_n(-200, jointScore) 

Top200_Pi = results %>%
  filter(SummaryStat == "piPerBin") %>%
  top_n(-200, jointScore)

####Generate high density regions and compute credible intervals and find mode 
#If you do not supply the density of the parameter this package will smooth the posterior using default density, otherwise kernel density estimation will be done with it's own algorithm where the default kernel bandwidth h is selected using the algorithm of Samworth and Wand (2010)
cols = c('SimNeThree','SimTimeTwo','SimNeTwo','SimTimeOne','SimNeOne')
credibleIntervals = data.frame()

for (Param in cols) {
  df = hdr(x = Top200[,Param], prob = c(95), den = density(Top200[,Param]))
  ParamCredInt = cbind.data.frame(Param, ceiling(df$mode), round(df$hdr, 3))
  #hdr.den(x = Top200[,Param], prob = c(95), den = density(Top200[,Param]))
  credibleIntervals = rbind.data.frame(credibleIntervals, ParamCredInt)
}

credibleIntervals$Param = c("Ancient population size", "Ancient bottleneck (generations ago)", "Intermidiate population size","Recent bottleneck (generations ago)", "Current population size") 
colnames(credibleIntervals) = c("Parameter", "Point Estimate", "Lower Bound (95% CI)", "Upper Bound (95% CI)")
ABCSummaryTable = ggtexttable(credibleIntervals, rows = NULL)


#Plotting Function for posterior
plotPosterior <- function(title, xAxisTitle, parameterPrior, colorPrior, parameterPosterior, colorPosterior, modePosterior ){
  ggplot() + 
    geom_density(data = sampPriors, aes(x = parameterPrior), colour=colorPrior,size=1, outline.type = "full") +
    geom_density(data = Top200, aes(x = parameterPosterior), colour=colorPosterior, outline.type = "full") +
    geom_vline(xintercept = modePosterior, colour="purple") +
    ggtitle(title) +
    xlab(xAxisTitle) +
    theme_bw() +
    theme(axis.text.x = element_text(size  = 20), 
          axis.text.y = element_text(size = 20), 
          axis.title=element_text(size=24, face = "bold"), 
          plot.title = element_text(size=24, hjust = 0.5, face = "bold"))
}

plotPosteriorNe1 = plotPosterior(paste("Mode of Posterior = ", credibleIntervals[5,2]), expression(paste(N[CURRENT])), sampPriors$Ne1, "blue", Top200$SimNeOne, "black", credibleIntervals[5,2])

plotPosteriorNe2 = plotPosterior(paste("Mode of Posterior = ", credibleIntervals[3,2]), expression(paste(N[INTERMIDIATE])), sampPriors$Ne2, "blue", Top200$SimNeTwo, "black", credibleIntervals[3,2])

plotPosteriorNe3 = plotPosterior(paste("Mode of Posterior = ", credibleIntervals[1,2]), expression(paste(N[ANCIENT])), sampPriors$Ne3, "blue", Top200$SimNeThree, "black", credibleIntervals[1,2])

plotPosteriorTbot1 = plotPosterior(paste("Mode of Posterior = ", credibleIntervals[4,2]), expression(paste(TBOT[RECENT])), sampPriors$Tbot, "blue", Top200$SimTimeOne, "black", credibleIntervals[4,2])

plotPosteriorTbot2 = plotPosterior(paste("Mode of Posterior = ", credibleIntervals[2,2]), expression(paste(TBOT[ANCIENT])), sampPriors$Tbot2, "blue", Top200$SimTimeTwo, "black", credibleIntervals[2,2])

plotPosteriorTogether = ggarrange(plotPosteriorNe1  + ylab(""),
                                  plotPosteriorTbot1  + ylab(""), 
                                  plotPosteriorNe2  + ylab(""), 
                                  plotPosteriorTbot2  +ylab(""), 
                                  plotPosteriorNe3  + ylab(""),
                                  nrow = 5, 
                                  ncol = 1,
                                  align = 'v')

#Parse simulated data
HetBins = Top200 %>%
  select(starts_with('V')) %>% 
  mutate_all(funs(sum(as.numeric(as.character(.)), na.rm = TRUE))) %>%
  distinct() %>%
  gather() %>%
  mutate(proportion = value / sum(value),
         bin = as.numeric(1:n()-1),
         Data = "Simulated") 

#Parse simulated data
piBins = Top200_Pi %>%
  select(starts_with('V')) %>% 
  mutate_all(funs(sum(as.numeric(as.character(.)), na.rm = TRUE))) %>%
  distinct() %>%
  gather() %>%
  mutate(proportion = value / sum(value),
         bin = as.numeric(1:n()-1),
         Data = "Simulated") 

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
  mutate(PIranges = as.numeric(as.character(cut(windowedPi, seq(from = 0, to = 0.025, by = 0.0005), include.lowest = T, labels = c(0:49))))) %>%
  group_by(PIranges) %>%
  count() %>%
  ungroup() %>%
  mutate(propSites = n/sum(n),
         Data = "Empirical") %>%
  dplyr::rename(bin = PIranges)

#Combine data for Seg Sites
subsetReal = HetsRealData %>% select("bin","propSites", "Data")
subsetSim = HetBins %>% select("bin","proportion", "Data")
matchCols = c("bin","propSites","Data")
colnames(subsetReal) <- matchCols
colnames(subsetSim) <- matchCols
comboDF = rbind.data.frame(subsetReal,subsetSim)

CompDataS = ggplot(data=comboDF %>% filter(bin <= 13), 
                   aes(x=bin, y=propSites,fill=Data)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() + 
  labs(x=expression(italic(S)), y="Proportion of data") +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=20)) +
  scale_fill_manual(values = c("Empirical"="blue","Simulated"="red")) +
  scale_x_continuous(limits=c(-1, 13), breaks=seq(from=0,to=12))

CompDataS_zoom = ggplot(data=comboDF %>% filter(bin > 0 & bin < 13), 
                        aes(x=bin, y=propSites,fill=Data)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x=expression(italic(S)), y="Proportion of data") +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) + 
  scale_fill_manual(values = c("Empirical"="blue","Simulated"="red")) +
  ylim(0,0.025) + 
  scale_x_continuous(limits=c(0.5, 13), breaks=seq(1:12)) 

#Combine data for Pi
subsetReal_pi = PiRealData %>% select("bin","propSites", "Data")
subsetSim_pi = piBins %>% select("bin","proportion", "Data")
matchCols = c("bin","propSites","Data")
colnames(subsetReal_pi) <- matchCols
colnames(subsetSim_pi) <- matchCols
comboDF_pi = rbind.data.frame(subsetReal,subsetSim)

#create data frame with ranges 
ranges = cut(realData$windowedPi, seq(from = 0, to = 0.025, by = 0.0005), include.lowest = T)
bin = as.numeric(0:49)
rangesDF = as.character(levels(ranges)) %>% 
  cbind(bin) %>% 
  as.data.frame()
names(rangesDF)[1] = "ranges"

#Add intervals to comboDF
comboDF_pi$intervals = rangesDF$ranges[match(comboDF_pi$bin,rangesDF$bin)]
comboDF_pi$intervals = factor(comboDF_pi$intervals, levels = rangesDF$ranges) #order the intervals

###Plot Pi
CompDataPi = ggplot(data=comboDF_pi %>% filter(bin <= 12), aes(x=intervals, y=propSites,fill=Data)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() + 
  labs(x=expression(pi), y="Proportion of data") +
  theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 0.5), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.title = element_blank(), 
        legend.text=element_text(size=20)) + 
  scale_fill_manual(values = c("Empirical"="blue","Simulated"="red"))

CompDataPi_zoom = ggplot(data=comboDF_pi %>% filter(bin > 0 & bin <= 12), aes(x=intervals, y=propSites,fill=Data)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() + 
  labs(x=expression(pi), y="Proportion of data") +
  theme(axis.text.x = element_text(size  = 20, angle = 90, hjust = 0.5), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size=24), 
        plot.title = element_text(size=24, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=20)) + 
  scale_fill_manual(values = c("Empirical"="blue","Simulated"="red")) +
  ylim(0,0.025)

#Plot Pi and S together
plotAll = ggarrange(CompDataS, 
                    CompDataS_zoom + ylab(""),
                    CompDataPi,
                    CompDataPi_zoom + ylab(""), 
                    align = 'v',
                    nrow = 2, 
                    ncol = 2,
                    common.legend = TRUE,
                    legend = "bottom")

ggarrange(plotAll, plotPosteriorTogether, ncol = 2)

###Plot joint distribution
meanSumStatEmpirical = realData %>%
  summarise_at(c("windowedPi","windowedS"), sum, na.rm=TRUE) %>%
  dplyr::rename(totalS = windowedS, totalPi = windowedPi) %>%
  as.data.frame()

results$Status = ifelse(results$jointScore %in% Top200$jointScore, "Accepted", "Rejected")

meanSumStatSimulated = read.table("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/results_3epoch/ResultsABC_3epoch_NAnc200k_200Ksims/allReplicates_meanSumStats_3epoch.txt", col.names = c("totalS","totalPi","Fname")) 


test = cbind.data.frame(meanSumStatSimulated[,c("totalS", "totalPi")], results)


DensityS = ggplot() + 
  geom_density(data = test, aes(x=totalS, colour=Status)) +
  scale_colour_manual(values=c(Accepted="blue", Rejected="gray30")) + 
  theme_bw() + 
  geom_vline(xintercept = as.numeric(meanSumStatEmpirical$totalS))

DensityPi = ggplot() +
  geom_density(data = test, aes(x=totalPi, colour=Status)) +
  scale_colour_manual(values=c(Accepted="blue", Rejected="gray30")) +
  theme_bw() +
  geom_vline(xintercept = as.numeric(meanSumStatEmpirical$totalPi))


VisualizeABC = ggplot() +
  geom_point(data = test %>% filter(Status == "Rejected" & totalS <=2500), 
             aes(x=totalS,
                 y=totalPi, 
                 colour = Status,
                 alpha = Status)) +
  scale_colour_manual(values=c(Accepted="blue", Rejected="gray60")) +
  scale_alpha_manual(values=c(Rejected = 0.2, Accepted = 1),guide=F) +
  geom_point(data = test %>% filter(Status == "Accepted" & totalS <= 2500), 
             aes(x=totalS,
                 y=totalPi, 
                 colour = Status,
                 alpha = Status)) +
  geom_point(data = meanSumStatEmpirical, 
             aes(x=totalS,
                 y=totalPi), 
             shape=8, 
             color="red", 
             size=2) + 
  theme_bw() +
  labs(x="S", y=expression(pi)) +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size=24, face = "bold"), 
        plot.title = element_text(size=24, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=20))



#Option 1 to generate an hpd set of level p, based on a sample x from the posterior
#hpd<-function(x,p){
#  dx<-density(x)
#  md<-dx$x[which.max(dx$y)]
#  px<-dx$y/sum(dx$y)
#  pxs<--sort(-px)
#  ct<-min(pxs[cumsum(pxs)< p])
#  list(hpdr=range(dx$x[px>=ct]),mode=md)
#}
#hpd(Top200$SimNeOne, 0.95)