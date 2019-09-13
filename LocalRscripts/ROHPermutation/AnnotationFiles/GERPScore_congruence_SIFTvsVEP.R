#load library 
library(tidyverse)

#load annotations
annotation = read_delim("~/Documents/DogProject_Clare/LocalRscripts/annotateVCF/AllChroms/annotatedGTwithVEP_AW10_allChroms.txt", delim = "\t")

#split into two data frames 
VEP = annotation %>%
  select(CHROM, POS, ANNOT, IMPACT, GERPScore) %>%
  na.omit()

#summarize the data
TP = VEP %>%
  filter(GERPScore > 4 & IMPACT == "HIGH" ) %>%
  count(IMPACT) %>%
  select(n) %>%
  as.numeric()
TN = VEP %>%
  filter(GERPScore < 2 & IMPACT == "LOW" ) %>%
  count(IMPACT)%>%
  select(n) %>%
  as.numeric()
FN = VEP %>%
  filter(GERPScore > 4 & IMPACT == "LOW" ) %>%
  count(IMPACT)%>%
  select(n) %>%
  as.numeric()
FP = VEP %>%
  filter(GERPScore < 2 & IMPACT == "HIGH" ) %>%
  count(IMPACT)%>%
  select(n) %>%
  as.numeric()
PPV = TP/(TP+FP)
SENSITIVITY = TP/(TP+FN)
SPECIFICITY = TN/(TN+FP)
sumStats_VEP = cbind.data.frame(FP,TN,TP,FN, PPV, SENSITIVITY, SPECIFICITY)

####NOW USE SIFT ANNOTATIONS####
SIFTscore = annotation %>%
  select(CHROM, POS, ANNOT, SIFT, GERPScore) %>%
  na.omit()

#summarize the data
TP = SIFTscore %>%
  filter(GERPScore > 4 & SIFT == "deleterious" ) %>%
  count(SIFT) %>%
  select(n) %>%
  as.numeric()
TN = SIFTscore %>%
  filter(GERPScore < 2 & SIFT == "benign" ) %>%
  count(SIFT)%>%
  select(n) %>%
  as.numeric()
FN = SIFTscore %>%
  filter(GERPScore > 4 & SIFT == "benign" ) %>%
  count(SIFT)%>%
  select(n) %>%
  as.numeric()
FP = SIFTscore %>%
  filter(GERPScore < 2 & SIFT == "deleterious" ) %>%
  count(SIFT)%>%
  select(n) %>%
  as.numeric()
PPV = TP/(TP+FP)
SENSITIVITY = TP/(TP+FN)
SPECIFICITY = TN/(TN+FP)
sumStats_SIFT = cbind.data.frame(FP,TN,TP,FN, PPV, SENSITIVITY, SPECIFICITY)

#Merge together
sumStats_VEP$Annot = "VEP"
sumStats_SIFT$Annot = "SIFT"

merged = rbind.data.frame(sumStats_VEP, sumStats_SIFT)

write.table(merged, file = "contigencyTableGERPScore_SIFTvsVEPannot.txt",col.names = T, row.names = F, quote = F, sep = "\t")

