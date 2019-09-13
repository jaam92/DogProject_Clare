#Set Libraries and working directory
library(ggplot2)
library(scales)
library(dplyr)
library(data.table)
library(mgsub)
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

botStart = ggplot(EW10_25iter,aes(x=generationsLeft,y=Ne, colour=INDV))+
  geom_step(stat="identity")+
  geom_vline(xintercept = 2400)+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(paste("EW10 PSMC Results\nmu= ",mu,"\ngeneration time = ",gen," yrs/gen",sep=""))+
  xlab("Generations Ago")+
  ylab("Ne")+
  scale_y_log10(labels=comma,breaks=c(1000,10000,100000,1000000))+
  scale_x_log10(breaks=c(100,1000,10000,100000),labels=comma)+
  theme(axis.text.x = element_text(size  = 20, vjust=1, hjust=1), axis.text.y = element_text(size  = 20), axis.title=element_text(size=20),legend.title=element_blank(), legend.text=element_text(size=18),legend.position="bottom")

#Plot results
ggplot(results25,aes(x=LeftYears,y=Ne, colour=INDV))+
  geom_step(stat="identity")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(paste("PSMC Results\nmu= ",mu,"\ngeneration time = ",gen," yrs/gen",sep=""))+
  xlab("Years Ago")+
  ylab("Ne")+
  scale_y_log10(labels=comma,breaks=c(1000,10000,100000,1000000))+
  scale_x_log10(breaks=c(100,1000,10000,100000,1000000),labels=comma)+
  theme(legend.position="bottom")

  
