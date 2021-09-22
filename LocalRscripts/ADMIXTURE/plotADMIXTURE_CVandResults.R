#Load library
library(ggplot2)
library(dplyr)
library(mgsub)
library(reshape2)

#Load files
CVError = read.delim("~/Documents/DogProject_Clare/LocalRscripts/ADMIXTURE/CrossValidationError.out")
QFile = read.table("~/Documents/DogProject_Clare/LocalRscripts/ADMIXTURE/prunedData_masterFile_allChroms.6.Q")
FamFile = read.table("~/Documents/DogProject_Clare/LocalRscripts/ADMIXTURE/prunedData_masterFile_allChroms.fam")
Individuals = read.table("~/Documents/DogProject_Clare/LocalRscripts/Dogs2Keep.txt")

#Plot the Error
ggplot(CVError, aes(x=K_subpops,y=CV_error)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() #low point is K = 6

#Make data frame of admixture for K=6
orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf", "Isle Royale")

plotResults = cbind.data.frame(FamFile$V1, QFile) %>%
  rename(INDV = 'FamFile$V1') %>%
  mutate(Population = substr(INDV,1,2), #make population column
         Population = mgsub(Population,pattern = c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement=c("Border Collie","Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf","Isle Royale")),
         Population = factor(Population, levels = orderPops),
         INDV = factor(INDV, levels = Individuals$V1)) %>% 
  melt(id=c("INDV","Population")) #melt data to plot

#plot results
(ADMIXTURE = ggplot(plotResults, aes(x=INDV,y=value,fill=variable)) + 
    geom_bar(stat="identity") + 
    theme_bw() + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x="Individual", y="Ancestry Proportion") + 
    theme(axis.text.x = element_text(size = 10, angle = 90, hjust=0.5, vjust = 0), 
          axis.text.y = element_text(size = 18), 
          plot.title=element_text(size=26, face = "bold", hjust=0.5), 
          axis.title=element_text(size=20), 
          legend.position = "none"))


#Average ancestry per cluster
AvgAncestry = cbind.data.frame(FamFile$V1, QFile) %>%
  rename(INDV = 'FamFile$V1') %>%
  mutate(Population = substr(INDV,1,2), #make population column
         Population = mgsub(Population,pattern = c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement=c("Border Collie","Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf","Isle Royale")),
         Population = factor(Population, levels = orderPops),
         INDV = factor(INDV, levels = Individuals$V1)) %>%
  group_by(Population) %>%
  select(-c(INDV)) %>%
  summarise_all(mean) 