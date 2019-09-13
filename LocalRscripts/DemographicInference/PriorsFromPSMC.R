#Set Libraries and working directory
library(ggplot2)
library(scales)
library(dplyr)
library(data.table)
library(mgsub)
library(cowplot)
setwd("~/Documents/DogProject_Clare/LocalRscripts/outputMSMC")
#Load vcftools outputs
fnames = list.files(pattern = "\\.final.txt$")
results = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE, idcol = "FileName")#Read all the files and create a FileName column to store filenames
mu=4.96e-9#Tanya Autosomal Mutation Rate for Neutral cm >0.4 and 3years/gen
gen=3 #(years/gen)
# to convert: from the msmc guide
# MSMC outputs times and rates scaled by the mutation rate per basepair per generation. First, scaled times are given in units of the per-generation mutation rate. This means that in order to convert scaled times to generations, divide them by the mutation rate. In humans, we used mu=1.25e-8 per basepair per generation.To convert generations into years, multiply by the generation time, for which we used 30 years.

# To get population sizes out of coalescence rates, first take the inverse of the coalescence rate, scaledPopSize = 1 / lambda. Then divide this scaled population size by 2*mu (yes, this factor 2 is different from the time scaling, sorry).
results$Ne <- (1/results$lambda)/(2*mu) # note the factor of 2! (not in time scaling) confirmed correct: https://github.com/stschiff/msmc-tools/blob/master/plot_utils.py
results$LeftYears <- gen*(results$left_time_boundary/mu)
results$RightYears <- gen*(results$right_time_boundary/mu)
results$INDV = mgsub(results$FileName, pattern = c(".*_", ".final.txt"),replacement = c("",""))
results$Population = substr(results$INDV,1,2)
results$numIterations = mgsub(results$FileName, pattern = c("_[^_]*$", "msmc2_output_iter"),replacement = c("","") )
results25 = results %>% filter(numIterations == "25")
#Get parameters for priors 
EW10_25iter = results25 %>% filter(INDV == "EW10") %>% mutate(generationsLeft = LeftYears/3)

modelPSMC = ggplot(EW10_25iter,aes(x=generationsLeft,y=Ne, colour=INDV))+
  geom_step(stat="identity")+
  geom_vline(xintercept = 2400)+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(paste("EW10 PSMC Results\nmu= ",mu,"\ngeneration time = ",gen," yrs/gen",sep=""))+
  xlab("Generations Ago")+
  ylab("Ne")+
  scale_y_log10(labels=comma,breaks=c(1000,10000,100000,1000000))+
  scale_x_log10(breaks=c(100,1000,10000,100000),labels=comma)+
  theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), axis.title=element_text(size=20),plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position=c(0.8, 0.8))

#Plot the priors
set.seed(1)
samps = seq(1:100000)
Ne1 = runif(samps, min = 900, max = 10000)
Ne2 = runif(samps, min = 11000, max = 40000)
Tbot = runif(samps, min = 1700, max = 2400)

plotNe1 = ggplot() + geom_density(aes(x=Ne1), colour="blue") + ggtitle(expression(paste(Ne[1]))) + xlab("unif(900,10000)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotNe2 = ggplot() + geom_density(aes(x=Ne2), colour="red") + ggtitle(expression(paste(Ne[2]))) + xlab("unif(11000,40000)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotTbot = ggplot() + geom_density(aes(x=Tbot), colour="green") + ggtitle(expression(paste(T[bottleneck]))) + xlab("unif(1700,2400)") + theme_bw()+ theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

plotPriorsTogether = plot_grid(plotNe1,plotNe2+ylab(""),plotTbot+ylab(""),nrow = 1, align = 'h')
plot_grid(modelPSMC, plotPriorsTogether, labels = "AUTO", ncol = 1, align = 'v', axis = 'r') # aligning vertically along the right axis