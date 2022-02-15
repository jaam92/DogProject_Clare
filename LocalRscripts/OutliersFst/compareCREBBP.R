#Load libraries and set working directory
library(ggrepel)
library(ggpubr)
setwd("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst")

#Load gene files
genes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1.bed", sep = "\t")
gene_names = read.table("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt")

#Load up functions and parse gene sets 
source("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FunctionsForCREBBPComparison.R")

#Generate Dataframes for pi
dfBC_pi = makeDataFrames("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/BC_allSites.sites.pi") %>%
  mutate(Population = "BC")
dfEW_pi = makeDataFrames("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/EW_allSites.sites.pi") %>%
  mutate(Population = "EW")
dfTM_pi = makeDataFrames("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/TM_allSites.sites.pi") %>%
  mutate(Population = "TM")
dfAW_pi = makeDataFrames("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/AW_allSites.sites.pi") %>%
  mutate(Population = "AW")

#Compute Pi across all genes
BC_piPerGene = computePI(dfBC_pi, "BC")
EW_piPerGene = computePI(dfEW_pi, "EW")
TM_piPerGene = computePI(dfTM_pi, "TM")
AW_piPerGene = computePI(dfAW_pi, "AW")

#Finescale evaluation of Pi for CREBBP
BC_CREBBP = computePI_CREBBP(dfBC_pi, "BC")
EW_CREBBP = computePI_CREBBP(dfEW_pi, "EW")
TM_CREBBP = computePI_CREBBP(dfTM_pi, "TM")
AW_CREBBP = computePI_CREBBP(dfAW_pi, "AW")

mergedDF_pi = rbind.data.frame(BC_piPerGene,TM_piPerGene,EW_piPerGene, AW_piPerGene) %>%
  mutate(chrom = gsub("chr", "", GeneSet$chrom[match(GeneName, GeneSet$AbbrevName)]),
         fullPopName = case_when(Population == "AW" ~ "Arctic Wolf",
                                 Population == "BC" ~ "Border Collie",
                                 Population == "TM" ~ "Tibetan Mastiff",
                                 Population == "EW" ~ "Ethiopian Wolf")) 

#data frame with CREBBP
CREBBP = mergedDF_pi %>%
  filter(GeneName == "CREBBP") 

allSitesDF_pi = rbind.data.frame(dfBC_pi, dfTM_pi, dfEW_pi, dfAW_pi) %>%
  filter(GeneName == "CREBBP") %>%
  mutate(fullPopName = case_when(Population == "AW" ~ "Arctic Wolf",
                                 Population == "BC" ~ "Border Collie",
                                 Population == "TM" ~ "Tibetan Mastiff",
                                 Population == "EW" ~ "Ethiopian Wolf"))


####Generate dataframes and compute Fst across all genes
BC_FstPerGene = compFST("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/EW_vs_BC_allSites_rmNAN.weir.fst") %>%
  mutate(Population = "BC")
TM_FstPerGene = compFST("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/EW_vs_TM_allSites_rmNAN.weir.fst") %>%
  mutate(Population = "TM")
AW_FstPerGene = compFST("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/EW_vs_AW_allSites_rmNAN.weir.fst") %>%
  mutate(Population = "AW")

mergedDF_Fst = rbind.data.frame(BC_FstPerGene,TM_FstPerGene,AW_FstPerGene) %>%
  mutate(chrom = gsub("chr", "", GeneSet$chrom[match(GeneName, GeneSet$AbbrevName)]),
         fullPopName = case_when(Population == "AW" ~ "Ethiopian wolf versus Arctic wolf",
                                 Population == "BC" ~ "Ethiopian wolf versus Border Collie",
                                 Population == "TM" ~ "Ethiopian wolf versus Tibetan Mastiff"))

#data frame with CREBBP Fst
CREBBP_Fst = mergedDF_Fst %>%
  filter(GeneName == "CREBBP") 

####Data frames for fixed sites sites can be fixed as follows:
  #ANC in comp group and DER in EW
  #ANC in EW and DER in comp group
CountPerGene_EWvAW = compFixedSites("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsAW_N6.txt") %>%
  filter(GeneName%in%EW_piPerGene$GeneName) %>% #only compare against gene sets used for computing pi
  mutate(Population = "AW")

CountPerGene_EWvTM = compFixedSites("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsTM_N6.txt") %>%
  filter(GeneName%in%EW_piPerGene$GeneName)  %>% 
  mutate(Population = "TM") %>%
  na.omit()

CountPerGene_EWvBC = compFixedSites("~/Documents/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsBC_N6.txt") %>%
  filter(GeneName%in%EW_piPerGene$GeneName) %>% 
  mutate(Population = "BC") %>%
  na.omit()

mergedDF_FixedSites = rbind.data.frame(CountPerGene_EWvAW, CountPerGene_EWvTM, CountPerGene_EWvBC) %>%
  mutate(fullPopName = case_when(Population == "AW" ~ "Ethiopian wolf versus Arctic wolf",
                        Population == "BC" ~ "Ethiopian wolf versus Border Collie",
                        Population == "TM" ~ "Ethiopian wolf versus Tibetan Mastiff")) %>%
  na.omit() 

#data frame with CREBBP fixed sites
CREBBP_FixedSites = mergedDF_FixedSites %>%
  filter(GeneName == "CREBBP")

####Start plotting
cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

####Plots for pi
dodge = position_dodge(width = 0.9)
ggplot(mergedDF_pi, aes(x=fullPopName,y=meanPI,colour=fullPopName)) + 
  geom_violin(size=1, position = dodge)  + 
  scale_colour_manual(name = "Population", values = cbPalette) + 
  ylab(expression("Mean" ~ pi ~ "per gene")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16))

#Histogram
ggplot() +
  geom_density(data = AW_piPerGene, aes(x=meanPI), colour ="gray25") +
  geom_density(data = TM_piPerGene, aes(x=meanPI), colour ="#CC79A7") +
  geom_density(data = BC_piPerGene, aes(x=meanPI), colour ="#009E73") +
  geom_density(data = EW_piPerGene, aes(x=meanPI), colour = "#D55E00") +
  xlab(expression("Mean" ~ pi ~ "per gene")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16))

#plot distribution with each species as facet
CREBBP$pvalues = c(paste0("p = ", round(digits = 3, x = nrow(BC_piPerGene[BC_piPerGene$meanPI<=0.0228, ])/nrow(BC_piPerGene)), sep=""), 
                   paste0("p = ", round(digits = 3, x = nrow(TM_piPerGene[TM_piPerGene$meanPI<=0.0417, ])/nrow(TM_piPerGene)), sep=""), 
                   paste0("p = ", round(digits = 3, x = nrow(EW_piPerGene[EW_piPerGene$meanPI<=0.00482, ])/nrow(EW_piPerGene)), sep=""), 
                   paste0("p = ", round(digits = 3, x = nrow(AW_piPerGene[AW_piPerGene$meanPI<=0.0174, ])/nrow(AW_piPerGene)), sep=""))
#uncomment these coordinates if only using counts
#CREBBP$x = c(0.2, 0.2, 0.1, 0.3) 
#CREBBP$y = c(290, 290, 890, 290)
CREBBP$x = c(0.25, 0.25, 0.25, 0.25) 
CREBBP$y = c(0.5, 0.5, 0.5, 0.5)

#plot distribution with each species as facet
ggplot(data = mergedDF_pi, aes(meanPI, group=fullPopName)) + 
  geom_histogram(aes(y = stat(count)/10045, fill=fullPopName), breaks = seq(0,1, by=0.05), color="black") + ####10045 is the number of genes in comparison
  scale_fill_manual(name = "Population", values = cbPalette) + 
  geom_vline(data = CREBBP, aes(xintercept = meanPI), colour="black", linetype = "dashed") + 
  geom_text(data = CREBBP, mapping = aes(x = x, y = y, label = pvalues), size = 6) +
  facet_wrap(~fullPopName) +
  labs(x=expression("Mean" ~ pi ~ "per gene"), y = "Proportion of genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.position = "none")


#CREBBP only plot as dots
ggplot(CREBBP, aes(x=fullPopName,y=meanPI,colour=fullPopName)) + 
  geom_point(size=1)  + 
  scale_colour_manual(name = "Population", values = cbPalette) + 
  ylab(expression("Mean" ~ pi ~ "CREBBP")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16))

#plots for windowed pi
CREBBP_1Mb_pi = piWindowedAboutGene(dfEW,"chr6", 36412372, 38531822, "1Mb", "CREBBP")
GRB2_1Mb_pi = piWindowedAboutGene(dfEW,"chr9", 4109244, 6181472, "1Mb", "GRB2")
MKL1_1Mb_pi = piWindowedAboutGene(dfEW,"chr10", 23599593, 25752076, "1Mb", "MKL1")
PYGB_1Mb_pi = piWindowedAboutGene(dfEW,"chr23", 426984, 2491240, "1Mb", "PYGB")
CDK8_1Mb_pi = piWindowedAboutGene(dfEW,"chr25", 12058947, 14178073, "1Mb", "CDK8")

#####Plotting for Fst
CREBBP_Fst$pvalues = c(paste0("p = ", round(digits = 3, x = nrow(BC_FstPerGene[BC_FstPerGene$FstNormSNPCount>=0.8919779, ])/nrow(BC_FstPerGene)), sep=""), 
                       paste0("p = ", round(digits = 3, x = nrow(TM_FstPerGene[TM_FstPerGene$FstNormSNPCount>=0.8616547, ])/nrow(TM_FstPerGene)), sep=""), 
                       paste0("p = ", round(digits = 3, x = nrow(AW_FstPerGene[AW_FstPerGene$FstNormSNPCount>=0.8865586, ])/nrow(AW_FstPerGene)), sep=""))
#CREBBP_Fst$x = c(0.25, 0.25, 0.2)
#CREBBP_Fst$y = c(240, 240, 240)
CREBBP_Fst$x = c(0.15, 0.15, 0.15)
CREBBP_Fst$y = c(0.12, 0.12, 0.12)

#plot distribution with each species as facet
#use sapply to get proportion per group
ggplot(mergedDF_Fst, aes(FstNormSNPCount, group=fullPopName)) + 
  geom_histogram(aes(y =stat(count)/sapply(PANEL, FUN=function(x) sum(count[PANEL == x])), 
                     fill=fullPopName), breaks = seq(0,1, by=0.05), color="black", fill="gray80") +
  geom_vline(data = CREBBP_Fst, aes(xintercept = FstNormSNPCount), colour="black", linetype="dashed") + 
  geom_text(data = CREBBP_Fst, mapping = aes(x = x, y = y, label = pvalues), size = 6) +
  facet_wrap(~fullPopName) +
  labs(x=expression("Mean" ~F[ST]~"per gene"), y="Frequency") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))

####Plots for derived allele counts
CREBBP_FixedSites$pvalues = c(paste0("p = ", round(digits = 3, x = nrow(CountPerGene_EWvAW[CountPerGene_EWvAW$FixedSites>=212, ])/nrow(CountPerGene_EWvAW)), sep=""), 
                              paste0("p = ", round(digits = 3, x = nrow(CountPerGene_EWvTM[CountPerGene_EWvTM$FixedSites>=121, ])/nrow(CountPerGene_EWvTM)), sep=""), 
                              paste0("p = ", round(digits = 3, x = nrow(CountPerGene_EWvBC[CountPerGene_EWvBC$FixedSites>=214, ])/nrow(CountPerGene_EWvBC)), sep=""))
#uncomment for counts
#CREBBP_FixedSites$x = c(600, 600, 600)
#CREBBP_FixedSites$y = c(1900, 1900, 1900)
CREBBP_FixedSites$x = c(600, 600, 600)
CREBBP_FixedSites$y = c(0.35, 0.35, 0.35)

#plot with populations as facets
ggplot(mergedDF_FixedSites, aes(FixedSites, group=fullPopName)) + 
  geom_histogram(aes(y=stat(count)/sapply(PANEL, FUN=function(x) sum(count[PANEL == x]))), breaks = seq(0,2000, by=50), color="black", fill="gray80") +
  geom_vline(data = CREBBP_FixedSites, aes(xintercept = FixedSites), colour="black", linetype="dashed") + 
  geom_text(data = CREBBP_FixedSites, mapping = aes(x = x, y = y, label = pvalues), size = 6) +
  facet_wrap(~fullPopName) +
  labs(x="Count of oppositely fixed sites (per gene)", y="Frequency") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))
