#Load Libraries
library(ggplot2)
library(dplyr)
library(IRanges)
library(GenomicRanges)
library(cowplot)

#Make all chroms a single data frame
setwd("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/InputData")
df = read.delim("EW4_allChroms_NeutralRegions_het_1000win_1000step.txt")

#Remove windows with no sites
#Window size is 1 million so we want coverage percentage based on 10^6
annotateDF = function(inputDF){
  inputDF %>%
    filter(sites_passing != 0) %>%
    mutate(chromo = as.numeric(as.character(gsub("chr", "", chromo))), 
           window_end = window_start + 1000, 
           TwentyPerCov = ifelse(sites_total >= 200, "1", "0"),
           FiftyPerCov = ifelse(sites_total >= 500, "1", "0"), 
           EightyPerCov = ifelse(sites_total >= 800, "1", "0")) 
}

Annotated = annotateDF(df) 


#Get ready to plot
# Prepare the dataset
ProcessDF = function(inputDF){ inputDF %>% 
    
    # Compute chromosome size
    group_by(chromo) %>% 
    summarise(chr_len=max(window_start)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(inputDF, ., by=c("chromo"="chromo")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chromo, window_start) %>%
    mutate( newWinStart=window_start+tot) }

#Preprocess and create raw data frame and
#Make new data frames based on coverage cutoffs``
plotEW = ProcessDF(Annotated) %>% 
  filter(TwentyPerCov == 1) #Coverage used to select prior was 20% callable sites

#Plot
#Generate x axis with one of the data frames 
axisdf = plotEW %>% 
  group_by(chromo) %>% 
  summarize(center=( max(newWinStart) + min(newWinStart) ) / 2 )

plotFunction = function(dataFrame, indiv) {
  indivHet = ggplot() + 
    geom_bar(data = dataFrame, aes(x=newWinStart, y=indiv/sites_total, color=as.factor(chromo)),stat = "identity", lwd=0.5) +
    scale_color_manual(values = rep(c("#E7298A", "#7570B3"), 38 )) +
    #custom X axis:
    scale_x_continuous(label = axisdf$chromo, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0)) + # remove space between plot area and x axis
    #custom theme
    labs(x = "chromosome", y = "heterozygosity(per bp)") +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none")
  return(indivHet)
}

(neutralhets_EW4 = plotFunction(plotEW, plotEW$hets_EW4))

#Plot distribution of number of hets 
(DensityPlotHets = ggplot() + 
    geom_density(data=plotEW, aes(x=hets_EW4), fill="blue") + 
    theme_bw() + 
    xlab("Number of hets per 1Kb") + 
    theme(axis.text.x = element_text(size  = 20), 
          axis.text.y = element_text(size = 20), 
          axis.title=element_text(size=20), 
          plot.title = element_text(size=20, hjust = 0.5),
          legend.title=element_blank(), 
          legend.text=element_text(size=18),legend.position="bottom"))

(DensityPlotHetsZoom = ggplot() + 
    geom_density(data=plotEW, aes(x=hets_EW4), fill="blue") +
    theme_bw() + 
    xlab("Number of hets per 1Kb") + 
    theme(axis.text.x = element_text(size  = 20), 
          axis.text.y = element_text(size = 20), 
          axis.title=element_text(size=20), 
          plot.title = element_text(size=20, hjust = 0.5),
          legend.title=element_blank(), 
          legend.text=element_text(size=18),legend.position="bottom"))

#Plot proportion with histogram
totalHets = dim(plotEW)[1]
histogramHets = plotEW %>%
  group_by(hets_EW4) %>%
  count() %>%
  mutate(propBin = n/totalHets) %>%
  ungroup()

(HistogramPlotHets = ggplot() + 
    geom_bar(data=histogramHets, aes(x=hets_EW4, y=propBin), stat = "identity") +
    theme_bw() + 
    labs(x="Count of Heterozygotes per 1Kb", y="Proportion of Total Hets") + 
    theme(axis.text.x = element_text(size  = 20), 
          axis.text.y = element_text(size = 20), 
          axis.title=element_text(size=20), 
          plot.title = element_text(size=20, hjust = 0.5),
          legend.title=element_blank(), 
          legend.text=element_text(size=18),legend.position="bottom"))

Heterozygosity = sum(plotEW$hets_EW4)/sum(plotEW$sites_total)
Heterozygosity
