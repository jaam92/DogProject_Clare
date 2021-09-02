#Load libraries and set working directory
library(ggrepel)
library(ggpubr)
setwd("~/DogProject_Clare/LocalRscripts/OutliersFst")

#Load gene files
genes = read.delim("~/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1.bed", sep = "\t")
gene_names = read.table("~/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt")

#Load up functions and parse gene sets 
source("~/DogProject_Clare/LocalRscripts/OutliersFst/FunctionsForCREBBPComparison.R")

#Generate Dataframes for pi
dfBC_pi = makeDataFrames("~/DogProject_Clare/LocalRscripts/OutliersFst/BC_allSites.sites.pi") %>%
  mutate(Population = "BC")
dfEW_pi = makeDataFrames("~/DogProject_Clare/LocalRscripts/OutliersFst/EW_allSites.sites.pi") %>%
  mutate(Population = "EW")
dfTM_pi = makeDataFrames("~/DogProject_Clare/LocalRscripts/OutliersFst/TM_allSites.sites.pi") %>%
  mutate(Population = "TM")
dfAW_pi = makeDataFrames("~/DogProject_Clare/LocalRscripts/OutliersFst/AW_allSites.sites.pi") %>%
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
BC_FstPerGene = compFST("~/DogProject_Clare/LocalRscripts/OutliersFst/EW_vs_BC_allSites_rmNAN.weir.fst") %>%
  mutate(Population = "BC")
TM_FstPerGene = compFST("~/DogProject_Clare/LocalRscripts/OutliersFst/EW_vs_TM_allSites_rmNAN.weir.fst") %>%
  mutate(Population = "TM")
AW_FstPerGene = compFST("~/DogProject_Clare/LocalRscripts/OutliersFst/EW_vs_AW_allSites_rmNAN.weir.fst") %>%
  mutate(Population = "AW")

mergedDF_Fst = rbind.data.frame(BC_FstPerGene,TM_FstPerGene,AW_FstPerGene) %>%
  mutate(chrom = gsub("chr", "", GeneSet$chrom[match(GeneName, GeneSet$AbbrevName)]),
         fullPopName = case_when(Population == "AW" ~ "Ethiopian Wolf versus Arctic Wolf",
                                 Population == "BC" ~ "Ethiopian Wolf versus Border Collie",
                                 Population == "TM" ~ "Ethiopian Wolf versus Tibetan Mastiff"))

#data frame with CREBBP Fst
CREBBP_Fst = mergedDF_Fst %>%
  filter(GeneName == "CREBBP") 

####Data frames for fixed sites
#Ethiopian wolf versus arctic
CountPerGene_EWvAW = compFixedSites("~/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsAW_N9.txt") %>%
  filter(GeneName%in%EW_piPerGene$GeneName) #only compare against gene sets used for computing pi

CountPerGene_EWvTM = compFixedSites("~/DogProject_Clare/LocalRscripts/OutliersFst/FixedSites_EWvsTM_N9.txt") %>%
  filter(GeneName%in%EW_piPerGene$GeneName) #only compare against gene sets used for computing pi


####Plots for pi
cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

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
ggplot() + 
  geom_histogram(data = mergedDF_pi, aes(x=meanPI), bins=100) +
  geom_vline(data = CREBBP, aes(xintercept = meanPI), colour="purple") + 
  facet_wrap(~fullPopName, scale = "free") +
  xlab(expression("Mean" ~ pi ~ "per gene")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))

#plot as dots
ggplot(CREBBP, aes(x=fullPopName,y=meanPI,colour=fullPopName)) + 
  geom_point() +
  scale_colour_manual(name = "Population", values = cbPalette) + 
  ylab(expression("Mean" ~ pi ~ "CREBBP")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16))

#Pi per site
ggplot(allSitesDF_pi, aes(x=POS, y=PI, colour=fullPopName)) +
  geom_point() + 
  scale_colour_manual(name = "Population", values = cbPalette) + 
  labs(x = "Position" ,y=expression(pi)) + 
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
#plot distribution with each species as facet
ggplot() + 
  geom_histogram(data = mergedDF_Fst, aes(x=FstNormSNPCount), bins=100) +
  geom_vline(data = CREBBP_Fst, aes(xintercept = FstNormSNPCount), colour="purple") + 
  facet_wrap(~fullPopName, scale = "free") +
  xlab(expression(F[ST])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))


####Plots for derived allele counts
ggplot(CountPerGene_EWvAW, aes(x=n)) + 
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 212, colour="purple") + #Count for CREBBP
  ggtitle(paste0("p-value = ", round(digits = 3, x = nrow(CountPerGene_EWvAW[CountPerGene_EWvAW$n>212, ])/nrow(CountPerGene_EWvAW)))) +
  labs(x="Number of Fixed Sites Per Gene (EW vs AW)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16))

#Ethiopoan wolf versus tibetan mastiff
ggplot(CountPerGene_EWvTM, aes(x=n)) + 
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 72, colour="purple") + #Count for CREBBP
  ggtitle(paste0("p-value = ", round(digits = 3, x = nrow(CountPerGene_EWvTM[CountPerGene_EWvTM$n>72, ])/nrow(CountPerGene_EWvTM)))) +
  labs(x="Number of Fixed Sites Per Gene (EW vs TM)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16))
