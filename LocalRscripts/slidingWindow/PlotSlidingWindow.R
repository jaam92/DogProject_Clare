#Load Libraries
library(tidyverse)
library(data.table)
library(IRanges)
library(GenomicRanges)
library(ggpubr)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/slidingWindow")

#Function will annotate files so that have coverage percentages 
#Window size is 1 million so we want coverage percentage based on 10^6
annotateDF = function(inputDF){
  inputDF %>%
    filter(sites_passing != 0) %>%
    mutate(chromo = as.numeric(as.character(gsub("chr", "", chromo))), 
           window_end = window_start + 10^6, 
           TwentyPerCov = ifelse(sites_total >= 2e+05, "1", "0"),
           FiftyPerCov = ifelse(sites_total >= 5e+05, "1", "0"), 
           SixtyPerCov = ifelse(sites_total >= 6e+05, "1", "0")) 
}

#This function will make a new data frame and adjust windows
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

#This function implements the previous two functions to create our final data frame for plotting

makeFinalFiles <- function(filesToJoin){
  Fnames = lapply(Sys.glob(filesToJoin), read.delim) #grab all chroms
  df = rbindlist(Fnames) #merge all chroms
  dfAnnot = annotateDF(df) %>% 
    filter(FiftyPerCov == 1) #use only windows with at least 50 percent coverage if you increase TM doesn't have very many windows
  processedDF = ProcessDF(dfAnnot)
  return(processedDF)
}

#Generate dataframes from all sliding windows
EW_Annotated = makeFinalFiles("TenEW_jointcalled_chr*_het_1000000win_1000000step.txt")
AW_Annotated = makeFinalFiles("All15AW_jointcalled_chr*_het_1000000win_1000000step.txt")
IR_Annotated = makeFinalFiles("TenIR_jointcalled_chr*_het_1000000win_1000000step.txt")
LB_Annotated = makeFinalFiles("TenLB_jointcalled_chr*_het_1000000win_1000000step.txt")
PG_Annotated = makeFinalFiles("All15PG_jointcalled_chr*_het_1000000win_1000000step.txt")
TM_Annotated = makeFinalFiles("TenTM_jointcalled_chr*_het_1000000win_1000000step.txt")
BC_Annotated = makeFinalFiles("TenBC_jointcalled_chr*_het_1000000win_1000000step.txt")

#Function to plot data for each individual
plotFunction = function(dataFrame, indiv, color1, color2) {
  #Generate x axis with any one of the data frames 
  axisdf = dataFrame %>% 
    group_by(chromo) %>% 
    summarize(center=( max(newWinStart) + min(newWinStart) ) / 2 )
  #Now plot with the axis
  indivHet = ggplot() + 
    geom_bar(data = dataFrame, aes(x=newWinStart, y=indiv/sites_total, color=as.factor(chromo)),stat = "identity", lwd=0.5) +
    scale_color_manual(values = rep(c(color1, color2), 38 )) +
    #custom X axis:
    scale_x_continuous(label = axisdf$chromo, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0)) + # remove space between plot area and x axis
    labs(x = "Chromosome", y = "Heterozygosity(per bp)") +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none")
  return(indivHet)
}

#Plot Ethiopian Wolves

#hets_EW1 = plotFunction(EW_Annotated, EW_Annotated$hets_EW1, "#D55E00", "bisque")
#hets_EW2 = plotFunction(EW_Annotated, EW_Annotated$hets_EW2, "#D55E00", "bisque")
hets_EW3 = plotFunction(EW_Annotated, EW_Annotated$hets_EW3, "#D55E00", "bisque")
#hets_EW4 = plotFunction(EW_Annotated, EW_Annotated$hets_EW4, "#D55E00", "bisque")
#hets_EW5 = plotFunction(EW_Annotated, EW_Annotated$hets_EW5, "#D55E00", "bisque")
#hets_EW6 = plotFunction(EW_Annotated, EW_Annotated$hets_EW6, "#D55E00", "bisque")
#hets_EW7 = plotFunction(EW_Annotated, EW_Annotated$hets_EW7, "#D55E00", "bisque")
#hets_EW8 = plotFunction(EW_Annotated, EW_Annotated$hets_EW8, "#D55E00", "bisque")
#hets_EW9 = plotFunction(EW_Annotated, EW_Annotated$hets_EW9, "#D55E00", "bisque")
#hets_EW10 = plotFunction(EW_Annotated, EW_Annotated$hets_EW10, "#D55E00", "bisque")

#ggarrange(hets_EW1,hets_EW2,hets_EW3,hets_EW4,hets_EW5,hets_EW5,hets_EW6,hets_EW8,hets_EW9,hets_EW10, nrow = 5, ncol = 2, align='hv')


#Plot IR 
#IR relatedness (1,6,9,4 are unknown anestry; 10 & 2 parent-offspring; 7 & 8 parent-offspring; 5 & 3 are inbred af )

#hets_IR1 = plotFunction(IR_Annotated, IR_Annotated$hets_IR1,"#FF9326", "steelblue")
#hets_IR2 = plotFunction(IR_Annotated, IR_Annotated$hets_IR2,"#FF9326", "steelblue")
hets_IR3 = plotFunction(IR_Annotated, IR_Annotated$hets_IR3,"#FF9326", "steelblue")
#hets_IR4 = plotFunction(IR_Annotated, IR_Annotated$hets_IR4,"#FF9326", "steelblue")
#hets_IR5 = plotFunction(IR_Annotated, IR_Annotated$hets_IR5,"#FF9326", "steelblue")
#hets_IR6 = plotFunction(IR_Annotated, IR_Annotated$hets_IR6,"#FF9326", "steelblue")
#hets_IR7 = plotFunction(IR_Annotated, IR_Annotated$hets_IR7,"#FF9326", "steelblue")
#hets_IR8 = plotFunction(IR_Annotated, IR_Annotated$hets_IR8,"#FF9326", "steelblue")
#hets_IR9 = plotFunction(IR_Annotated, IR_Annotated$hets_IR9,"#FF9326", "steelblue")
#hets_IR10 = plotFunction(IR_Annotated, IR_Annotated$hets_IR10,"#FF9326", "steelblue")

#ggarange(hets_IR1,hets_IR6,hets_IR9,hets_IR4,hets_IR10,hets_IR7,hets_IR2,hets_IR8,hets_IR5,hets_IR3, nrow = 5, ncol = 2, align='hv')

#Plot one example from the others
hets_AW13 = plotFunction(AW_Annotated, AW_Annotated$hets_AW13,"gray25", "cyan3")
hets_PG3 = plotFunction(PG_Annotated, PG_Annotated$hets_PG3,"mediumpurple4", "aquamarine2")
hets_BC3 = plotFunction(BC_Annotated, BC_Annotated$hets_BC3,"#009E73", "mistyrose")
hets_LB3 = plotFunction(LB_Annotated, LB_Annotated$hets_LB3,"gold3", "firebrick2")
hets_TM3 = plotFunction(TM_Annotated, TM_Annotated$hets_TM3,"#CC79A7", "lightblue2")