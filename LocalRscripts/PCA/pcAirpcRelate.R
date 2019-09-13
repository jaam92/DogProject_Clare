#Load necessary packages
library(tidyverse)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(ggrepel)
library(mgsub)

#Set working directory and load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/PCA")
#snpgdsBED2GDS(bed.fn = "masterFile_allChroms_updateFID.bed", bim.fn = "masterFile_allChroms_updateFID.bim", fam.fn = "masterFile_allChroms_updateFID.fam", out.gdsfn = "masterFile_allChroms_updateFID.gds")

#Open files and pre
gds = snpgdsOpen("masterFile_allChroms_updateFID.gds")
sampIds = read.gdsn(index.gdsn(gds, "sample.id")) #grab sample ids 
famIds = substr(sampIds,1,2) #make family ids

#LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
snpset = snpgdsLDpruning(gds, sample.id=sampIds, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
pruned = unlist(snpset, use.names=FALSE)

#Run KING
KING = snpgdsIBDKING(gds, sample.id = sampIds, snp.id=pruned, autosome.only = F)

#Make GRM
KINGmat = KING$kinship %>%
  as.matrix()  
colnames(KINGmat) = KING$sample.id 
rownames(KINGmat) = KING$sample.id 

#Extract pairwise kinship estimates and IBS0
dfIBS = snpgdsIBDSelection(KING) %>%
  mutate(pair = paste(ID1, ID2, sep = "_"))

snpgdsClose(gds)

#Plot kinship data
#ggplot(dfIBS, aes(IBS0, kinship)) +
#  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
#  geom_point(alpha=0.5) +
#  ylab("kinship estimate") +
#  theme_bw()

#ggplot(dfIBS, aes(IBS0, kinship)) +
#  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
#  geom_point(alpha=0.5) +
#  ylab("kinship estimate") +
#  ylim(-0.5,0.5) +
#  theme_bw() #zoom in on relateds

#ggplot(dfIBS, aes(IBS0, kinship, label = pair)) +
#  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
#  geom_point(alpha=0.5) +
#  ylab("kinship estimate") +
#  ylim(-0.5,0.5) +
#  theme_bw() +
#  geom_text_repel(data = subset(dfIBS, kinship > 0.20),
#                  nudge_y = 0.5 - subset(dfIBS, kinship > 0.20)$kinship,
#                  segment.size  = 0.2,
#                  segment.color = "grey50",
#                  direction = "x") #label the highly related inidividuals

#Open GDS data 
genoFile = GdsGenotypeReader(filename = "masterFile_allChroms_updateFID.gds")#read in GDS data
genoData = GenotypeData(genoFile)#create a GenotypeData class object

##Partition data into relateds and unrelated (less than first cousins)
set.seed(100)
sampset = pcairPartition(kinobj = KINGmat, kin.thresh=2^(-9/2), divobj = KINGmat, div.thresh=2^(-9/2)) 
unrelsID = sampset$unrels

#Run PC-AiR
canidpcAir = pcair(genoData, kinobj = KINGmat, kin.thresh=2^(-9/2), divobj = KINGmat, div.thresh=2^(-9/2), unrel.set = unrelsID, snp.include = pruned, autosome.only = F)
summary(canidpcAir)


#Prep GDS for PC-Relate
genoData = GenotypeBlockIterator(genoData, snpBlock = 20000)#take input from PC-AiR and convert to Genotype block iterator

#Run PC-Relate
canidpcRelate = pcrelate(genoData, pcs = canidpcAir$vectors[,1:2], training.set = canidpcAir$unrels, ibd.probs = FALSE) #use first 2 pcs to correct kinship for population structure (aka ancestry)
pcRelateMat = pcrelateToMatrix(canidpcRelate, scaleKin = 1) #convert pcrelate output to GRM and don't scale kinship 

#Redo PCA with ancestry adjusted kinship 
correctedCanidpcAir = pcair(genoData, kinobj= pcRelateMat, kin.thresh=2^(-9/2), divobj= KINGmat, div.thresh=2^(-9/2), snp.include = pruned, autosome.only = F) #use ancestry adjusted pcrelate GRM in the next pca
summary(correctedCanidpcAir)

snpgdsClose(gds)

#Plot ancestry adjusted PCA
pcs = correctedCanidpcAir$vectors
pc.df = as.data.frame(pcs)
names(pc.df) = paste0("PC", 1:ncol(pcs))
pc.df$sample.id = row.names(pcs)
pc.df$Population = substr(pc.df$sample.id, 1, 2)
pc.df$Population = mgsub(pc.df$Population, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))

cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

PC1vPC2 = ggplot(pc.df, aes(x=PC1, y=PC2, color=Population)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Population") + 
  theme_bw() + 
  labs(y=bquote('PC2' ~'('~.(round(correctedCanidpcAir$values[2], digits = 3))~'%'~')'), x=bquote('PC1'~'('~.(round(correctedCanidpcAir$values[1], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

PC2vPC3 = ggplot(pc.df, aes(x=PC2, y=PC3, color=Population)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Population") + 
  theme_bw() + 
  labs(y=bquote('PC3' ~'('~.(round(correctedCanidpcAir$values[3], digits = 3))~'%'~')'), x=bquote('PC2'~'('~.(round(correctedCanidpcAir$values[2], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))