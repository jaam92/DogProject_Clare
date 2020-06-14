#Load Libraries
source("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/R_rainclouds.R")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgsub)

#Plotting fxn and color palette
cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

plotFxn = function(colOfInterest, axisTitle) {
  RaincloudWithBoxPlot = ggplot(PlotDF, aes(x=Population, y=colOfInterest, colour=Population)) +
    geom_flat_violin(size=1, position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE) +
    geom_point(aes(x = as.numeric(Population)-.15, y = colOfInterest, colour = Population),position = position_jitter(width = .05), size = 1, shape = 20) +
    geom_boxplot(aes(as.numeric(Population), y = colOfInterest, fill = Population), outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
    coord_flip() +
    guides(fill = FALSE, colour = FALSE) +
    scale_colour_manual(values = cbPalette) +
    scale_fill_manual(values = cbPalette) + 
    labs(x="Population",y=paste0(axisTitle)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size=20, angle = 45, vjust=0.7), 
          axis.text.y = element_text(size=20), 
          plot.title=element_text(size=24, hjust = 0.5, face = "bold"), 
          axis.title=element_text(size=24),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=20))
  return(RaincloudWithBoxPlot)
}

#Read file in 
df = read.delim("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AllChroms/GTAnnotationCountResults_May2019_DogProjClare.txt", stringsAsFactors = F)

####Make Plotting Data Frame
PlotDF = df %>% 
  mutate(CallableSites = LineCount - Missing)

PlotDF$Population = mgsub(PlotDF$Population, 
                          pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))
orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale")
PlotDF$Population = factor(PlotDF$Population, levels = orderPops)

####Make Everything proportional 
PropPlotDF = PlotDF[,c(1:11, 20:25)] %>%
  mutate_at(vars(LOF_CountAlleles: PutNeu_CountVariants), list(~./PlotDF$CallableSites)) #Just run replace on plotting fxns with regular df to make everything proportional

####Make Scaled Count Alleles
scaleCalls = mean(PlotDF$CallableSites)
ScaledPlotDF = PropPlotDF %>%
  mutate_at(vars(LOF_CountAlleles: PutNeu_CountVariants), list(~.*scaleCalls))

pSY = plotFxn(ScaledPlotDF$SY_CountAlleles, "Synonymous Derived Alleles") 
pNS = plotFxn(ScaledPlotDF$NS_CountAlleles, "Nonsynonymous Derived Alleles") 
pLOF = plotFxn(ScaledPlotDF$LOF_CountAlleles, "Loss of function Derived Alleles") 

###Plot Supplementary Figures putativeley neutral and deleterious
CountDerHom_PutDel = plotFxn(PlotDF$PutDel_CountDerHom, "Count Derived Homozygotes") 
CountDerHom_PutNeu = plotFxn(PlotDF$PutNeu_CountDerHom, "Count Derived Homozygotes")

CountVar_PutDel = plotFxn(PlotDF$PutDel_CountVariants, "Count Variants") 
CountVar_PutNeu = plotFxn(PlotDF$PutNeu_CountVariants, "Count Variants")

CountAllele_PutDel = plotFxn(PlotDF$PutDel_CountAlleles, "Count Alleles") 
CountAllele_PutNeu = plotFxn(PlotDF$PutNeu_CountAlleles, "Count Alleles")

#Neutral
putNeu = ggarrange( CountDerHom_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Homozygotes"),
                    CountVar_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Variants"),
                    CountAllele_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Alleles"),
                    align = 'hv',
                    labels = c("A", "B", "C"),
                    hjust = -1,
                    nrow = 1)
putNeuAnnot = annotate_figure(putNeu, left = text_grob("Neutral", size=24, face="bold",rot = 90, hjust = 0.5))

#Deleterious
putDel = ggarrange( CountDerHom_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    CountVar_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    CountAllele_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    align = 'hv',
                    labels = c("D", "E", "F"),
                    hjust = -1,
                    nrow = 1)
putDelAnnot = annotate_figure(putDel, left = text_grob("Deleterious", size=24, face="bold",rot = 90, hjust = 0.5))

#plot deleterious and neutral
ggarrange(putNeuAnnot, putDelAnnot, nrow = 2)

#arrange the three scaled derived allele count plots in a single row
ggarrange( pSY,
           pNS + theme(axis.title.y=element_blank()),
           pLOF + theme(axis.title.y=element_blank()),
           align = 'hv',
           labels = c("A", "B", "C"),
           hjust = -1,
           nrow = 1)
