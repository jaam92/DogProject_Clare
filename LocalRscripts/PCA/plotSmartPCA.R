#Load necessary packages
library(tidyverse)
library(mgsub)

#Set working directory and load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/PCA")

#Read files in
pc.df = read_delim("~/Documents/DogProject_Clare/LocalRscripts/PCA/prunedData_masterFile_allChroms_smartPCA.snpeigs", delim = "\t") %>%
  mutate(Population = substr(Individual, 1, 2))
pc.variance = read_table("~/Documents/DogProject_Clare/LocalRscripts/PCA/prunedData_masterFile_allChroms_smartPCA.eigs", col_names = c("Variance"))

#Reformat Population  Column
pc.df$Population = mgsub(pc.df$Population, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))

#Plot PCA
cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

PC1vPC2 = ggplot(pc.df, aes(x=PC1, y=PC2, color=Population)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Population") + 
  theme_bw() + 
  labs(y=bquote('PC2' ~'('~.(round(as.numeric(pc.variance[2,1]), digits = 3))~'%'~')'), x=bquote('PC1'~'('~.(round(as.numeric(pc.variance[1,1]), digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

PC2vPC3 = ggplot(pc.df, aes(x=PC2, y=PC3, color=Population)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Population") + 
  theme_bw() + 
  labs(y=bquote('PC3' ~'('~.(round(as.numeric(pc.variance[3,1]), digits = 3))~'%'~')'), x=bquote('PC2'~'('~.(round(as.numeric(pc.variance[2,1]), digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
