#Load Libraries
source("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/R_rainclouds.R")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgsub)

#Plotting fxn and color palette
cbPalette = c("Arctic wolf" = "gray25", "Ethiopian wolf" = "#D55E00",  "Isle Royale wolf" = "steelblue", "border collie" = "#009E73", "labrador retriever" = "gold3", "pug" = "mediumpurple4", "tibetan mastiff" = "#CC79A7")

plotFxn = function(dataFrame, colOfInterest, axisTitle) {
  RaincloudWithBoxPlot = ggplot(dataFrame, aes(x=Population, y=colOfInterest, colour=Population)) +
    geom_flat_violin(size=1, position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE) +
    geom_point(aes(x = as.numeric(Population)-.15, y = colOfInterest, colour = Population),position = position_jitter(width = .05), size = 1, shape = 20) +
    geom_boxplot(aes(as.numeric(Population), y = colOfInterest, fill = Population), outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
    coord_flip() +
    guides(fill = "none", colour = "none") +
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
df = read.delim("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AllChroms/GTAnnotationCountResults_June2020_DogProjClare.txt", stringsAsFactors = F)

####Make Plotting Data Frame
PlotDF = df %>% 
  mutate(CallableSites = LineCount - Missing)

PlotDF$Population = mgsub(PlotDF$Population, 
                          pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("border collie", "labrador retriever", "pug", "tibetan mastiff", "Arctic wolf",  "Ethiopian wolf", "Isle Royale wolf"))
orderPops = c("border collie", "labrador retriever", "pug", "tibetan mastiff", "Arctic wolf",  "Ethiopian wolf", "Isle Royale wolf")
PlotDF$Population = factor(PlotDF$Population, levels = orderPops)

####Make Everything proportional 
PropPlotDF = PlotDF[,c(1:11, 20:31)] %>%
  mutate_at(vars(LOF_CountAlleles: PutNeu_CountVariants), list(~./PlotDF$CallableSites)) #Just run replace on plotting fxns with regular df to make everything proportional

####Make Scaled Count Alleles
scaleCalls = mean(PlotDF$CallableSites)
ScaledPlotDF = PropPlotDF %>%
  mutate_at(vars(LOF_CountAlleles: PutNeu_CountVariants), list(~.*scaleCalls))

###Plot Supplementary Figures putativeley neutral and deleterious
CountDerHom_PutDel = plotFxn(PlotDF, PlotDF$PutDel_CountDerHom, "Count Derived Homozygotes") 
CountDerHom_PutNeu = plotFxn(PlotDF, PlotDF$PutNeu_CountDerHom, "Count Derived Homozygotes")

CountVar_PutDel = plotFxn(PlotDF, PlotDF$PutDel_CountVariants, "Count Variants") 
CountVar_PutNeu = plotFxn(PlotDF, PlotDF$PutNeu_CountVariants, "Count Variants")

CountAllele_PutDel = plotFxn(PlotDF, PlotDF$PutDel_CountAlleles, "Count Alleles") 
CountAllele_PutNeu = plotFxn(PlotDF, PlotDF$PutNeu_CountAlleles, "Count Alleles")

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
CountDerHom_NS = plotFxn(ScaledPlotDF, ScaledPlotDF$NS_CountDerHom, "Count Derived Homozygotes") 
CountDerHom_SY = plotFxn(ScaledPlotDF, ScaledPlotDF$SY_CountDerHom, "Count Derived Homozygotes")

CountVar_NS = plotFxn(ScaledPlotDF, ScaledPlotDF$NS_CountVariants, "Count Variants") 
CountVar_SY = plotFxn(ScaledPlotDF, ScaledPlotDF$SY_CountVariants, "Count Variants")

CountAllele_NS = plotFxn(ScaledPlotDF, ScaledPlotDF$NS_CountAlleles, "Count Alleles") 
CountAllele_SY = plotFxn(ScaledPlotDF, ScaledPlotDF$SY_CountAlleles, "Count Alleles")

#Neutral
SY = ggarrange( CountDerHom_SY + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Homozygotes"),
                    CountVar_SY + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Variants"),
                    CountAllele_SY + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Alleles"),
                    align = 'hv',
                    labels = c("A", "B", "C"),
                    hjust = -1,
                    nrow = 1)
SYAnnot = annotate_figure(SY, left = text_grob("Synonymous", size=24, face="bold",rot = 90, hjust = 0.5))

#Deleterious
NS = ggarrange( CountDerHom_NS + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    CountVar_NS + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    CountAllele_NS + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    align = 'hv',
                    labels = c("D", "E", "F"),
                    hjust = -1,
                    nrow = 1)
NSAnnot = annotate_figure(NS, left = text_grob("Nonsynonymous", size=24, face="bold",rot = 90, hjust = 0.5))

#plot deleterious and neutral
ggarrange(SYAnnot, NSAnnot, nrow = 2)



#Compute pvalues
PlotDF$Population = mgsub(as.character(PlotDF$Population), 
                          pattern =c("border collie", "labrador retriever", "pug", "tibetan mastiff", "Arctic wolf",  "Ethiopian wolf", "Isle Royale wolf"),
                          replacement =c("Dog", "Dog", "Dog", "Dog", "wolf", "EW", "wolf"))

means = PlotDF[,c(1:11, 20:31)] %>%
  group_by(Population) %>%
  summarise_at(vars(PutDel_CountAlleles:PutNeu_CountVariants), list(~mean(.x)))

#Deleterious Count Hom
793.3000/506.1220 #dog
793.3000/448.5417 #wolf
pairwise.wilcox.test(PlotDF$PutDel_CountDerHom, PlotDF$Population, p.adj = "bonf")$p.value

#Deleterious Count Vars
976.3000/861.2195 #dog
976.3000/976.0417 #wolf
pairwise.wilcox.test(PlotDF$PutDel_CountVariants, PlotDF$Population, p.adj = "bonf")$p.value


#Deleterious Count Alleles
1769.600/1367.341 #dog
1769.600/1424.583 #wolf
pairwise.wilcox.test(PlotDF$PutDel_CountAlleles, PlotDF$Population, p.adj = "bonf")$p.value
  