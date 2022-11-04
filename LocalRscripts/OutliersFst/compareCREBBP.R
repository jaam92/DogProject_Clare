#Load Libraries
library(tidyverse)
library(data.table)
library(ggrepel)
library(ggpubr)
source("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FunctionsForCREBBPComparison.R") #Load up functions and parse gene sets need to source this after reading in genes file

#Load gene files
genes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1.bed", sep = "\t")
gene_names = read.table("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt")

####Create data frames with pi
setwd("~/Documents/DogProject_Clare/LocalRscripts/ComputeandPlotPi/N6_Curr/perPopulationperGene")
fnames = list.files(pattern = "\\PerPopulation.pi$")


df = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE) %>%
  dplyr::rename(c("meanPi" = "pi", "GeneName"="gene"))

dfEW_pi = df %>%
  filter(population == "EW") %>%
  mutate(percentile = percent_rank(meanPi))

dfAW_pi = df %>%
  filter(population == "AW") %>%
  mutate(percentile = percent_rank(meanPi))

dfTM_pi = df %>%
  filter(population == "TM") %>%
  mutate(percentile = percent_rank(meanPi))

dfBC_pi = df %>%
  filter(population == "BC") %>%
  mutate(percentile = percent_rank(meanPi))


mergedDF_pi = rbind.data.frame(dfBC_pi, dfTM_pi, dfEW_pi, dfAW_pi) %>%
  mutate(fullPopName = case_when(population == "AW" ~ "Arctic wolf",
                                 population == "BC" ~ "border collie",
                                 population == "TM" ~ "tibetan mastiff",
                                 population == "EW" ~ "Ethiopian wolf"))

CREBBP = rbind.data.frame(dfBC_pi, dfTM_pi, dfEW_pi, dfAW_pi) %>%
  filter(GeneName == "CREBBP") %>%
  mutate(fullPopName = case_when(population == "AW" ~ "Arctic wolf",
                                 population == "BC" ~ "border collie",
                                 population == "TM" ~ "tibetan mastiff",
                                 population == "EW" ~ "Ethiopian wolf"))

#computeP-values
CREBBP$pvalues = c(paste0("p = ", round(digits = 3, x = nrow(dfBC_pi[dfBC_pi$meanPi<=4.016414e-05, ])/nrow(dfBC_pi)), sep=""), 
                   paste0("p = ", round(digits = 3, x = nrow(dfTM_pi[dfTM_pi$meanPi<=1.077712e-04, ])/nrow(dfTM_pi)), sep=""), 
                   paste0("p = ", round(digits = 3, x = nrow(dfEW_pi[dfEW_pi$meanPi<=4.265475e-05, ])/nrow(dfEW_pi)), sep=""), 
                   paste0("p = ", round(digits = 3, x = nrow(dfAW_pi[dfAW_pi$meanPi<=1.188775e-04, ])/nrow(dfAW_pi)), sep=""))

CREBBP$x = c(0.015, 0.015, 0.015, 0.015) 
CREBBP$y = c(0.25, 0.25, 0.25, 0.25)

####Create data frames for Fst 
setwd("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/") #same directory for fst and fixed sites
fnames = list.files(pattern = "\\mean.fst$")


mergedDF_fst = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FILENAME") %>%
  dplyr::rename(c("meanFst" = "fst", "GeneName"="gene", "Comparison"="FILENAME")) %>%
  mutate(Comparison = gsub("_allSites_rmNAN.weir.perGene.mean.fst", "", Comparison),
         fullPopName = case_when(Comparison == "EW_vs_AW" ~ "Ethiopian wolf vs\n Arctic wolf",
                                 Comparison == "EW_vs_BC" ~ "Ethiopian wolf vs\n border collie",
                                 Comparison == "EW_vs_TM" ~ "Ethiopian wolf vs\n tibetan mastiff",
                                 Comparison == "BC_vs_AW" ~ "border collie vs\n Arctic wolf",
                                 Comparison == "BC_vs_TM" ~ "border collie vs\n tibetan mastiff"))

split_mergedDF_fst = split(mergedDF_fst, f = mergedDF_fst$Comparison) 

CREBBP_fst = mergedDF_fst %>%
  filter(GeneName == "CREBBP") 

#compute p-values
CREBBP_fst$pvalues = c(paste0("p = ", round(digits = 3, x = nrow(split_mergedDF_fst$BC_vs_AW[split_mergedDF_fst$BC_vs_AW$meanFst>=0.07775059, ])/nrow(split_mergedDF_fst$BC_vs_AW)), sep=""), 
                       paste0("p = ", round(digits = 3, x = nrow(split_mergedDF_fst$BC_vs_TM[split_mergedDF_fst$BC_vs_TM$meanFst>=0.02993633, ])/nrow(split_mergedDF_fst$BC_vs_TM)), sep=""), 
                       paste0("p = ", round(digits = 3, x = nrow(split_mergedDF_fst$EW_vs_AW[split_mergedDF_fst$EW_vs_AW$meanFst>=0.88625552, ])/nrow(split_mergedDF_fst$EW_vs_AW)), sep=""), 
                       paste0("p = ", round(digits = 3, x = nrow(split_mergedDF_fst$EW_vs_BC[split_mergedDF_fst$EW_vs_BC$meanFst>=0.89182346, ])/nrow(split_mergedDF_fst$EW_vs_BC)), sep=""), 
                       paste0("p = ", round(digits = 3, x = nrow(split_mergedDF_fst$EW_vs_TM[split_mergedDF_fst$EW_vs_TM$meanFst>=0.86156290, ])/nrow(split_mergedDF_fst$EW_vs_TM)), sep=""))

CREBBP_fst$x = c(0.15, 0.15, 0.15, 0.15, 0.15)
CREBBP_fst$y = c(0.10, 0.10, 0.10, 0.10, 0.10)

####Data frames for fixed sites sites can be fixed as follows:
  #ANC in comp group and DER in EW
  #ANC in EW and DER in comp group
CountPerGene_EWvAW = compFixedSites("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsAW_N6.txt") %>%
  filter(GeneName%in%mergedDF_fst$GeneName) %>% #only compare against gene sets used for computing pi
  mutate(Population = "AW",
         percentile = percent_rank(FixedSites)) %>%
  na.omit()

CountPerGene_EWvTM = compFixedSites("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsTM_N6.txt") %>%
  filter(GeneName%in%mergedDF_fst$GeneName)  %>% 
  mutate(Population = "TM",
         percentile = percent_rank(FixedSites)) %>%
  na.omit()

CountPerGene_EWvBC = compFixedSites("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsBC_N6.txt") %>%
  filter(GeneName%in%mergedDF_pi$GeneName) %>% 
  mutate(Population = "BC",
         percentile = percent_rank(FixedSites)) %>%
  na.omit()

mergedDF_FixedSites = rbind.data.frame(CountPerGene_EWvAW, CountPerGene_EWvTM, CountPerGene_EWvBC) %>%
  mutate(fullPopName = case_when(Population == "AW" ~ "Ethiopian wolf vs\n Arctic wolf",
                        Population == "BC" ~ "Ethiopian wolf vs\n border collie",
                        Population == "TM" ~ "Ethiopian wolf vs\n tibetan mastiff")) %>%
  na.omit() 

###Fixed derived sites in ethiopian wolves only
fixedEW = read_delim("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FixedDerivedSites_EW_N9.txt", delim = "\t")
CREBBP_1Mb_fixedSites = countFixDerHomGene(fixedEW,"chr6", 36412372, 38531822, "1Mb", "CREBBP", "EW", "18")
MKL1_1Mb_fixedSites = countFixDerHomGene(fixedEW,"chr10", 23599593, 25752076, "1Mb", "MKL1", "EW", "18")
PYGB_1Mb_fixedSites = countFixDerHomGene(fixedEW,"chr23", 426984, 2491240, "1Mb", "PYGB", "EW", "18")
CDK8_1Mb_fixedSites = countFixDerHomGene(fixedEW,"chr25", 12058947, 14178073, "1Mb", "CDK8", "EW", "18")
ggarrange(CREBBP_1Mb_fixedSites, MKL1_1Mb_fixedSites, PYGB_1Mb_fixedSites, CDK8_1Mb_fixedSites, ncol = 1, nrow = 4, common.legend = TRUE, legend = "bottom")

#data frame with CREBBP fixed sites
CREBBP_FixedSites = mergedDF_FixedSites %>%
  filter(GeneName == "CREBBP")

####Start plotting
cbPalette = c("Arctic wolf" = "gray25", "Ethiopian wolf" = "#D55E00",  "Isle Royale wolf" = "steelblue", "border collie" = "#009E73", "labrador retriever" = "gold3", "pug" = "mediumpurple4", "tibetan mastiff" = "#CC79A7")
scaleFUN <- function(x) sprintf("%.1f", x) #round x and y axis digits 
dodge = position_dodge(width = 0.9)

###Plot distribution of pi with each population as facet
PiComps = ggplot(data = mergedDF_pi, aes(meanPi, group=fullPopName)) + 
  geom_histogram(aes(y = stat(count)/14931, fill=fullPopName), binwidth = 0.0005, closed="right", boundary=-0.5, color="black") + ####14931 is the number of genes in comparison
  scale_fill_manual(name = "Population", values = cbPalette) + 
  geom_vline(data = CREBBP, aes(xintercept = meanPi), colour="black", linetype="dashed", size = 2) + 
  geom_text(data = CREBBP, mapping = aes(x = x, y = y, label = pvalues), size = 16) +
  facet_wrap(~fullPopName, nrow = 4) +
  labs(x=expression("Mean" ~ pi ~ "(per gene)"), y = "Frequency") + 
  theme_bw()+ 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 42),
        strip.text = element_text(size = 42),
        panel.spacing.y = unit(10, "mm"),
        legend.position = "none")

###Plots for Fst
#plot distribution with each species as facet
FSTComps = ggplot(mergedDF_fst %>% filter(str_detect(Comparison, 'EW_vs_')), aes(meanFst, group=fullPopName)) + 
  geom_histogram(aes(y =stat(count)/9850), binwidth = 0.025, closed="right", boundary=-0.5, color="black", fill="gray80") +
  geom_vline(data = CREBBP_fst %>% filter(str_detect(Comparison, 'EW_vs_')), aes(xintercept = meanFst), colour="black", linetype="dashed", size = 2) + 
  geom_text(data = CREBBP_fst %>% filter(str_detect(Comparison, 'EW_vs_')), mapping = aes(x = x, y = y, label = pvalues), size = 16) +
  facet_wrap(~fullPopName, nrow = 3) +
  labs(x=expression("Mean" ~F[ST]~"(per gene)"), y="Frequency") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 42),
        strip.text = element_text(size = 42),
        panel.spacing.y = unit(10, "mm"))

###Plots for derived allele counts
CREBBP_FixedSites$pvalues = c(paste0("p = ", round(digits = 3, x = nrow(CountPerGene_EWvAW[CountPerGene_EWvAW$FixedSites>=212, ])/nrow(CountPerGene_EWvAW)), sep=""), 
                              paste0("p = ", round(digits = 3, x = nrow(CountPerGene_EWvTM[CountPerGene_EWvTM$FixedSites>=121, ])/nrow(CountPerGene_EWvTM)), sep=""), 
                              paste0("p = ", round(digits = 3, x = nrow(CountPerGene_EWvBC[CountPerGene_EWvBC$FixedSites>=214, ])/nrow(CountPerGene_EWvBC)), sep=""))
CREBBP_FixedSites$x = c(600, 600, 600)
CREBBP_FixedSites$y = c(0.35, 0.35, 0.35)

#plot with populations as facets
OppFixedSitesComp = ggplot(mergedDF_FixedSites, aes(FixedSites, group=fullPopName)) + 
  geom_histogram(aes(y=stat(count)/sapply(PANEL, FUN=function(x) sum(count[PANEL == x]))), breaks = seq(0, 2000, by=50), closed="right", boundary=-0.5, color="black", fill="gray80") +
  geom_vline(data = CREBBP_FixedSites, aes(xintercept = FixedSites), colour="black", linetype="dashed", size = 2) + 
  geom_text(data = CREBBP_FixedSites, mapping = aes(x = x, y = y, label = pvalues), size = 16) +
  facet_wrap(~fullPopName, nrow = 3) +
  labs(x="Count of oppositely fixed sites (per gene)", y="Frequency") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 42),
        strip.text = element_text(size = 42),
        panel.spacing.y = unit(10, "mm"))


x = ggplot() + theme_void() #blank plot between

pdf(file = "~/Desktop/Figure4.pdf", width = 80, height = 40)
print(ggarrange(x, x, FSTComps, PiComps, OppFixedSitesComp, CREBBP_1Mb_fixedSites, nrow = 2, ncol = 3))
dev.off()

#####Find pi and fst outliers that are as extreme as CREBBP
FSTandPiOutliers = mergedDF_fst %>% 
  filter(str_detect(Comparison, 'EW_vs_'), #only EW comparisons
         percentile >= 0.95) %>%
  group_by(GeneName) %>%
  count() %>%
  ungroup() %>%
  filter(n>2) %>%
  mutate(n = NULL,
         pi = dfEW_pi$meanPi[match(GeneName, dfEW_pi$GeneName)],
         FST_percentile_BC = split_mergedDF_fst$EW_vs_BC$percentile[match(GeneName, split_mergedDF_fst$EW_vs_BC$GeneName)]*100,
         FixedSites_percentile_BC = CountPerGene_EWvBC$percentile[match(GeneName, CountPerGene_EWvBC$GeneName)]*100,
         FST_percentile_AW = split_mergedDF_fst$EW_vs_AW$percentile[match(GeneName, split_mergedDF_fst$EW_vs_AW$GeneName)]*100,
         FixedSites_percentile_AW = CountPerGene_EWvAW$percentile[match(GeneName, CountPerGene_EWvAW$GeneName)]*100,
         FST_percentile_TM = split_mergedDF_fst$EW_vs_TM$percentile[match(GeneName, split_mergedDF_fst$EW_vs_TM$GeneName)]*100,
         FixedSites_percentile_TM = CountPerGene_EWvTM$percentile[match(GeneName, CountPerGene_EWvTM$GeneName)]*100) %>%
  filter(pi < 4.265475e-05) 

#write.table(FSTandPiOutliers, "~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FSTOutlierPiAtLeastAsExtAsCREBBP_SummaryStats_updatedPi.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

####Generate supplementary figure 5
#Comparison with EW
BC_vs_EW = mergedDF_fst %>%
  filter(str_detect(Comparison, 'EW_vs_BC')) %>%
  mutate(Comparison = "BC_vs_EW",
         fullPopName = "border collie vs\n Ethiopian wolf")

CREBBP_fst = CREBBP_fst %>%
  filter(Comparison == "EW_vs_BC") %>%
  mutate(Comparison = "BC_vs_EW",
         fullPopName = "border collie vs\n Ethiopian wolf") %>%
  rbind.data.frame(CREBBP_fst)
  
  

supp5 = rbind.data.frame(mergedDF_fst, BC_vs_EW) %>%
  mutate(Label = ifelse(percentile >= 0.95, "Top 5% Outliers","Remainder of Genes")) %>%
  na.omit()



EWComps_FST = ggplot(supp5 %>% filter(str_detect(Comparison, 'EW_vs')), 
                     aes(x=meanFst,y=numSNPs,colour=Label)) + 
  geom_point() +  
  scale_colour_manual(name = "Status", values=c("blue","red")) + 
  geom_vline(data = CREBBP_fst %>%
               filter(str_detect(Comparison, 'EW_vs')), 
             aes(xintercept = meanFst), colour="black", linetype="dashed", size = 1) +
  facet_wrap(~fullPopName, nrow = 4) +
  labs(x=expression("Mean" ~ F[ST] ~ "(per gene)"), y = "SNP Count") +  
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 42),
        strip.text = element_text(size = 42),
        legend.title=element_text(size=40), 
        legend.text=element_text(size=38),
        panel.spacing.y = unit(10, "mm"))



BCComps_FST = ggplot(supp5 %>% filter(str_detect(Comparison, 'BC_vs')), 
                     aes(x=meanFst,y=numSNPs,colour=Label)) + 
  geom_point() +  
  scale_colour_manual(name = "Status", values=c("blue","red")) + 
  geom_vline(data = CREBBP_fst %>%
               filter(str_detect(Comparison, 'BC_vs')), 
             aes(xintercept = meanFst), colour="black", linetype="dashed", size = 1) +
  facet_wrap(~fullPopName, nrow = 4) +
  labs(x=expression("Mean" ~ F[ST] ~ "(per gene)"), y = "SNP Count") +  
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 42),
        strip.text = element_text(size = 42),
        legend.title=element_text(size=40), 
        legend.text=element_text(size=38),
        panel.spacing.y = unit(10, "mm"))

ggarrange(EWComps_FST, BCComps_FST, ncol = 2, common.legend = TRUE, legend = "top")
