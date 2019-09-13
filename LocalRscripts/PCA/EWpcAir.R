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

#Open files and prep
gds = snpgdsOpen("masterFile_allChroms_updateFID.gds")
sampIds = read.gdsn(index.gdsn(gds, "sample.id")) #grab sample ids 
famIds = substr(sampIds,1,2) #make family ids
EWonly = sampIds[22:31] #pull EW samps

#LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
snpset = snpgdsLDpruning(gds, sample.id=EWonly, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
pruned = unlist(snpset, use.names=FALSE)

#Run KING
KING = snpgdsIBDKING(gds, sample.id = EWonly, snp.id=pruned, autosome.only = F)

#Make GRM
KINGmat = KING$kinship %>%
  as.matrix()  
colnames(KINGmat) = KING$sample.id 
rownames(KINGmat) = KING$sample.id 
snpgdsClose(gds)

#Reopen GDS data 
genoFile = GdsGenotypeReader(filename = "masterFile_allChroms_updateFID.gds")#read in GDS data
genoData = GenotypeData(genoFile)#create a GenotypeData class object

#Pull unrelateds (at most 1st cousins)
sampset = pcairPartition(kinobj = KINGmat, kin.thresh=2^(-9/2), divobj = KINGmat, div.thresh=2^(-9/2)) 
unrelsID = sampset$unrels

#Run pc air
EWpcAir = pcair(genoData, kinobj= KINGmat, kin.thresh=2^(-9/2), divobj= KINGmat, div.thresh=2^(-9/2), snp.include = pruned, unrel.set = unrelsID, autosome.only = F) 

snpgdsClose(gds)

#Plot ancestry adjusted PCA
pcs = EWpcAir$vectors
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
  labs(y=bquote('PC2' ~'('~.(round(EWpcAir$values[2], digits = 3))~'%'~')'), x=bquote('PC1'~'('~.(round(EWpcAir$values[1], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
