###Plot Figures

##import all the code we need for Figure 1
##Figure 1 uses ggpubr 
##dimensions for pdf are w=40 and h=12
library(ggpubr)
source("~/Documents/DogProject_Clare/LocalRscripts/PCA/IBS_Clustering.R") #Tree based on IBS
source("~/Documents/DogProject_Clare/LocalRscripts/ComputeandPlotPi/N6_Curr/plotPi_nonROHvsGenomeWide.R") #whole genome pi
source("~/Documents/DogProject_Clare/LocalRscripts/SFS/PlotSFS.R") #whole genome sfs

#Legend bottom page
sumstats = ggarrange(plotGW_newLegend, 
                     WolfOnlySFS + ylab("Proportion") + ggtitle(NULL), 
                     ncol = 1, 
                     nrow = 2, 
                     common.legend = T, 
                     legend = "right",  
                     labels = c("B", "C"), 
                     font.label = list(size = 24))
ggarrange(IBSTree, 
          sumstats, 
          ncol = 2, 
          nrow = 1, 
          align = 'hv', 
          labels = c("A"), 
          font.label = list(size = 24))

#Restart R session
#.rs.restartR()


##import all the code we need for Figure 2
##Figure 2 uses ggarrange
##dimensions for pdf are w=50 and h=22
library(ggpubr)
source("~/Documents/DogProject_Clare/LocalRscripts/slidingWindow/PlotSlidingWindow.R")
source("~/Documents/DogProject_Clare/LocalRscripts/ROH/ROH_Qual.R")
pairSlidingWindowandROH = function(slidingWinPlot, ROHPlot){
  ggarrange(slidingWinPlot + 
              ylim(0,0.005) + 
              ylab("") + 
              theme(axis.title.x=element_blank(),
                    axis.text = element_text(size = 16)),
            ROHPlot + 
              ylab("") +
              theme(axis.title.y=element_blank(), 
                    legend.position = "none",  
                    axis.text.x = element_text(size = 18),
                    axis.text.y = element_text(size = 16),
                    panel.border = element_blank(), 
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank()),
            ncol = 2,
            nrow = 1,
            align = 'hv')
}

BC = pairSlidingWindowandROH(hets_BC3, ROHRangePlot_BC)
LB = pairSlidingWindowandROH(hets_LB3, ROHRangePlot_LB)
PG = pairSlidingWindowandROH(hets_PG3, ROHRangePlot_PG)
TM = pairSlidingWindowandROH(hets_TM3, ROHRangePlot_TM)
AW = pairSlidingWindowandROH(hets_AW13, ROHRangePlot_AW)
EW = pairSlidingWindowandROH(hets_EW3, ROHRangePlot_EW)
IR = ggarrange(hets_IR3 + 
                 ylim(0,0.005) + 
                 ylab("") + 
                 theme(axis.title.x=element_text(size=20, face = "bold"),
                       axis.text = element_text(size = 16),
                       legend.position = "none"),
               ROHRangePlot_IR + 
                 xlab("") +
                 theme(axis.title.y=element_text(size=20, face = "bold"), 
                       axis.text.x = element_text(size = 18),
                       axis.text.y = element_text(size = 16),
                       panel.border = element_blank(), 
                       panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       legend.position = "none"),
               ncol = 2,
               nrow = 1,
               align = 'hv')

p2 = ggarrange(BC, LB, PG, TM, AW, EW, IR, nrow = 7)
annotP2 = annotate_figure(p2,
                          left = text_grob("Heterozygosity(per bp)", 
                                           size=20, 
                                           face = "bold",
                                           rot = 90))
cowplot::plot_grid(annotP2, legend, rel_widths = c(2,.1))

#Restart R session
#.rs.restartR()

##import all the code we need for Figure 3
##Figure 3 uses ggarrange
##dimensions for pdf are w=18 and h=20
library(ggpubr)
source("~/Documents/DogProject_Clare/LocalRscripts/DemographicInference/plotPosteriorMultSamp3epoch_NeCurr10Ne70KTbot1000compRealvsSim_10KbWindow.R")
ggarrange(ABCandDemog, plotPosteriorTogether)


#Figure 4
##Figure 4 uses ggarrange
##dimensions for pdf are w=36 and h=12
library(ggpubr)
source("~/Documents/DogProject_Clare/LocalRscripts/PRDM9/plotNJTree.R")
source("~/Documents/DogProject_Clare/LocalRscripts/Fst/PlotPairwise_EWvsAW.R")
source("~/Documents/DogProject_Clare/LocalRscripts/Fst/PlotPairwise_EWvsTM.R")

geneTrees = ggarrange(PRDM9Tree, 
                      GAPDHTree, 
                      align = 'hv', 
                      nrow = 2,
                      labels=c("A","B"), 
                      font.label = list(size = 24))

Fst = ggarrange(correlationFstSnpCount_EPAS1 + labs(x="",y=""), 
                correlationFstSnpCount_CREBBP + labs(x="",y=""), 
                nrow = 2, 
                ncol = 1, 
                align = 'hv', 
                labels=c("C","D"), 
                font.label = list(size = 24), 
                common.legend = TRUE, 
                legend = "right")

Fst_addAxes = annotate_figure(Fst, 
                  left = text_grob("Count SNPs", 
                                   size = 20, 
                                   face = "bold", 
                                   rot = 90),
                  bottom = text_grob(expression(F[ST]), 
                                   size = 20, 
                                   face = "bold"))

ggarrange(geneTrees, Fst_addAxes, ncol = 2)
