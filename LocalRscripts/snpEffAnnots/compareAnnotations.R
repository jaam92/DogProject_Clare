#Load libraries
library(ggplot2)
library(dplyr)
library(data.table)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/snpEffAnnots")
all_annots = read.table("~/Documents/DogProject_Clare/LocalRscripts/snpEffAnnots/anns_snpEff.txt")
canonical_annots = read.table("~/Documents/DogProject_Clare/LocalRscripts/snpEffAnnots/anns_snpEff_cannonicalOnly.txt")

#Merge and compare
compare = merge(all_annots, canonical_annots, by=c("V1","V2"), all = T)
compare$sameAnnots = ifelse(compare$V5.x == compare$V5.y, "1", "0")
table(compare$sameAnnots)

#check out annots
compare %>%
  filter(sameAnnots == 1) %>%
  select(V1, V2, V3.x, V4.x, V5.x) %>%
  count(V5.x)

differences = compare %>%
  filter(sameAnnots == 0) %>%
  as.data.frame()
  