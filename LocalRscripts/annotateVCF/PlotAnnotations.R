#Load Libraries
source("~/Documents/DogProject_Jaz/OMIA/R_rainclouds.R")
library(ggplot2)
library(cowplot)
library(dplyr)
library(mgsub)

#Read file in 
df = read.delim("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AllChroms/GTAnnotationCountResults_Dec2018_DogProjClare.txt")

#Plot Counts SYN Derived Alleles
PlotDF = df %>% 
  mutate(CallableSites = LineCount - Missing, 
         PropSY_CountAlleles = SY_CountAlleles/CallableSites,
         PropNS_CountAlleles = NS_CountAlleles/CallableSites,
         PropLOF_CountAlleles = LOF_CountAlleles/CallableSites)

scaleCalls = mean(PlotDF$CallableSites)


PlotDF$Population = mgsub_dict(PlotDF$Population,conversions =  list("BC" = "Border Collie", "LB" = "Labrador Retriever", "PG" = "Pug", "TM" = "Tibetan Mastiff", "AW" = "Arctic Wolf", "EW" = "Ethiopian Wolf", "IR" = "Isle Royale"))
orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale")
PlotDF$Population = factor(PlotDF$Population, levels = orderPops)

cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

(SYAlleles = ggplot(PlotDF, aes(x=Population, y=PropSY_CountAlleles*scaleCalls,colour=Population, label=ID)) + geom_boxplot()+ scale_colour_manual(name="Group", values = cbPalette) + labs(x="Group",y="Synonymous Derived Alleles") + theme_bw() + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)))


(NSAlleles = ggplot(PlotDF, aes(x=Population, y=PropNS_CountAlleles*scaleCalls,colour=Population, label=ID)) + geom_boxplot()+ scale_colour_manual(name="Group", values = cbPalette) + labs(x="Group",y="Nonsynonymous Derived Alleles") + theme_bw() + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)))


(LOFAlleles = ggplot(PlotDF, aes(x=Population, y=PropLOF_CountAlleles*scaleCalls,colour=Population, label=ID)) + geom_boxplot()+ scale_colour_manual(name="Group", values = cbPalette) + labs(x="Group",y="LOF Derived Alleles") + theme_bw() + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)))



pSY =ggplot(PlotDF,aes(x=Population, y=PropSY_CountAlleles*scaleCalls,colour=Population)) +
  geom_flat_violin(size=1, position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE) +
  geom_point(aes(x = as.numeric(Population)-.15, y = PropSY_CountAlleles*scaleCalls, colour = Population),position = position_jitter(width = .05), size = 1, shape = 20) +
  geom_boxplot(aes(as.numeric(Population), y = PropSY_CountAlleles*scaleCalls, fill = Population),outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  coord_flip() +
  guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) + 
  labs(x="Group",y="Synonymous Derived Alleles") + 
  theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24))


pNS = ggplot(PlotDF,aes(x=Population, y=PropNS_CountAlleles*scaleCalls,colour=Population)) +
  geom_flat_violin(size=1, position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE) +
  geom_point(aes(x = as.numeric(Population)-.15, y = PropNS_CountAlleles*scaleCalls, colour = Population),position = position_jitter(width = .05), size = 1, shape = 20) +
  geom_boxplot(aes(as.numeric(Population), y = PropNS_CountAlleles*scaleCalls, fill = Population),outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  coord_flip() +
  guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) + 
  labs(x="Group",y="Nonsynonymous Derived Alleles") +
  theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24))

pLOF =ggplot(PlotDF,aes(x=Population, y=PropLOF_CountAlleles*scaleCalls,colour=Population)) +
  geom_flat_violin(size=1, position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE) +
  geom_point(aes(x = as.numeric(Population)-.15, y = PropLOF_CountAlleles*scaleCalls, colour = Population),position = position_jitter(width = .05), size = 1, shape = 20) +
  geom_boxplot(aes(as.numeric(Population), y = PropLOF_CountAlleles*scaleCalls, fill = Population),outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  coord_flip() +
  guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) + 
  labs(x="Group",y="Loss of Function Derived Alleles") +
  theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24))


# arrange the three plots in a single row
plot_grid( pSY,
           pNS + theme(axis.title.y=element_blank()),
           pLOF + theme(axis.title.y=element_blank()),
           align = 'vh',
           labels = c("A", "B", "C"),
           hjust = -1,
           nrow = 1)
