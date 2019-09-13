#Load libraries
library(cowplot)
library(ggplot2)
library(mgsub)

#Load Files
setwd("~/Documents/DogProject_Clare/LocalRscripts/ROHPermutation/usingNSvsSY_Dec2018/")

#Plot Results Permutations
ROH_10Kb = read.csv("ratios_pvalues_pergroup_10KbROH.csv", stringsAsFactors = F)
ROH_1Mb = read.csv("ratios_pvalues_pergroup_1MbROH.csv", stringsAsFactors = F)
ROH_btwn10Kb1Mb = read.csv("ratios_pvalues_pergroup_btwnTenKbandOneMb.csv", stringsAsFactors = F) 

#prep data for plotting
ROH_10Kb$group = mgsub(ROH_10Kb$group, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))

ROH_1Mb$group = mgsub(ROH_1Mb$group, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))

ROH_btwn10Kb1Mb$group = mgsub(ROH_btwn10Kb1Mb$group, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))

orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale")

ROH_10Kb$group = factor(ROH_10Kb$group, levels = orderPops)
ROH_1Mb$group = factor(ROH_1Mb$group, levels = orderPops)
ROH_btwn10Kb1Mb$group = factor(ROH_btwn10Kb1Mb$group, levels = orderPops)

tenKb = ggplot(ROH_10Kb, aes(x=group, y=ratio, colour=type, group=group)) + 
  geom_line(aes(x=group, y=ratio), colour="black") + 
  geom_point(aes(x=group, y=ratio,size=pvalue<0.05)) + 
  coord_flip() + 
  theme_bw() + 
  geom_hline(yintercept = 1, colour="blue", linetype="dashed") + 
  scale_colour_manual(name = "Counting Method", values = c("alleleCopies"="#E69F00", "homozygotes"="#56B4E9","variants"="#009E73")) + 
  labs(x= "Group", y = "Odds-Ratio") + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))

oneMb = ggplot(ROH_1Mb, aes(x=group, y=ratio, colour=type, group=group)) + 
  geom_line(aes(x=group, y=ratio), colour="black") + 
  geom_point(aes(x=group, y=ratio,size=pvalue<0.05)) + 
  coord_flip() + 
  theme_bw() + 
  geom_hline(yintercept = 1, colour="blue", linetype="dashed") + 
  scale_colour_manual(name = "Counting Method", values = c("alleleCopies"="#E69F00", "homozygotes"="#56B4E9","variants"="#009E73")) + 
  labs(x= "Group", y = "Odds-Ratio") + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))

btwnTenKbandOneMb = ggplot(ROH_btwn10Kb1Mb, aes(x=group, y=ratio, colour=type, group=group)) + 
  geom_line(aes(x=group, y=ratio), colour="black") + 
  geom_point(aes(x=group, y=ratio,size=pvalue<0.05)) + 
  coord_flip() + 
  theme_bw() + 
  geom_hline(yintercept = 1, colour="blue", linetype="dashed") + 
  scale_colour_manual(name = "Counting Method", values = c("alleleCopies"="#E69F00", "homozygotes"="#56B4E9","variants"="#009E73")) + 
  labs(x= "Group", y = "Odds-Ratio") + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))

#Plot it
legend =  get_legend(tenKb + theme(legend.position="right"))
threePlots = plot_grid(tenKb + theme(legend.position="none"),
                      btwnTenKbandOneMb + theme(legend.position="none"),
                      oneMb + theme(legend.position="none"),
                      align = 'vh',
                      labels = c("A", "B","C"),
                      hjust = -1,
                      nrow = 1)
plot_grid(threePlots, legend, rel_widths = c(3, .3))
