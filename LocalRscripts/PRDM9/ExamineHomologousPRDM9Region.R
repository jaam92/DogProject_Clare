#Load Libraries
library(dplyr)
library(data.table)
library(mgsub)
library(reshape2)
library(tidyr)


#Load file and subset to homologous PRDM9 region
annots = read.delim("allIndivs_Chr5_annotatedGTwithVEP.txt") %>%
  filter(POS >= 63559618 & POS <= 63626468) %>%
  mutate(INDV = gsub(".*annotatedGTwithVEP_\\s*|_Chr.*", "", ID, perl = T),
         POPULATION = substr(INDV,1,2)) %>%
  select(INDV, POPULATION, POS, GT, ANNOT, IMPACT, withinROH_10Kb, withinROH_1Mb) 

#Look at annotations
IMPACT = annots %>%
  group_by(POPULATION) %>%
  count(IMPACT) %>%
  filter(IMPACT != "LOW") %>%
  arrange(desc(n))

ANNOT = annots %>%
  group_by(POPULATION) %>%
  count(ANNOT) %>%
  filter(ANNOT == "NS") %>%
  arrange(desc(n))

MISSING = annots %>%
  group_by(POPULATION) %>%
  count(GT) %>%
  filter(GT == "Missing")

TOTAL_SITES = annots %>%
  count(POPULATION)

#Put it together
TOTAL_SITES %>%
  rename("TOTAL_SITES" = n) %>%
  mutate(MISSING = MISSING$n[match(POPULATION, MISSING$POPULATION)],
         NONMISSING = if_else(is.na(MISSING), TOTAL_SITES, TOTAL_SITES - MISSING),
         COUNT_NS = ANNOT$n[match(POPULATION, ANNOT$POPULATION)],
         PROP_NS = (COUNT_NS/NONMISSING)*100) %>%
  arrange(desc(PROP_NS))

