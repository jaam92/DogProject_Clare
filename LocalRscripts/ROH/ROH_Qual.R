#Load libraries
library(data.table)
library(tidyverse)
library(hablar)
library(mgsub)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/ROH")
Individuals = read.table("~/Documents/DogProject_Clare/LocalRscripts/Dogs2Keep.txt")

#Load vcftools outputs
fnames = list.files(pattern = "\\ROHQC.LROH$")
df = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE, idcol = "FileName") #Read all the files and create a FileName column to store filenames

#Remove header lines and convert character columns to numeric
nonQCROH_allCanids = subset(df, df$CHROM != "CHROM") %>% 
  retype() %>% 
  mutate(Population = substr(INDV, 1, 2)) %>% 
  as.data.frame()

rm(df)#delete old df

#Plot the distribution of coverage of ROHs pre-filtering
preFilter = ggplot() + 
  geom_density(data = nonQCROH_allCanids, aes(x=Prop_Covered)) + 
  labs(x="Proportion of ROH Covered") +
  ggtitle("Pre-Filtering") + 
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24))
#print(preFilter)

#QC ROH Segments

#Old Filter 
##Will remove ROH if the proportion of ROH segment covered by snps is not within a SD of the mean
#z = data.table(nonQCROH_allCanids)
#z[,ToKeep := abs(nonQCROH_allCanids$Prop_Cov - mean(nonQCROH_allCanids$Prop_Cov)) < sd(nonQCROH_allCanids$Prop_Cov)][ToKeep  == TRUE] #create variable that identifies what to drop or keep 
#table(z$ToKeep, z$Population)/length(z$ToKeep) #Proportion that fall within drop or keep for each Population
#allCanidsTrueROH = subset(z, z$ToKeep == "TRUE")

#Filter chosen for project (ROH is at least 10kb)
allCanidsTrueROH = nonQCROH_allCanids %>%
  filter(Prop_Covered >= 0.4 & ROH_length != 1 & ROH_length != 0 & ROH_length >= 10000) %>% 
  select(-c(FileName)) %>%
  mutate(CHROM = paste("chr",CHROM,sep = ""))

#write.table(allCanidsTrueROH, file="FinalROHQCd_min10Kb_allIndivs_Dec2018_DogProjClare.LROH", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Plot the distribution of coverage of ROHs post-filtering
postFilter = ggplot() + 
  geom_density(data = allCanidsTrueROH, aes(x=Prop_Covered)) + 
  labs(x="Proportion of ROH Covered") + 
  ggtitle("Post-Filtering") + 
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24))
#print(postFilter)

#remove anything less than 1Mb so that FROH and ROH plots are the same ROH
#aggregate by individual and population and find FROH
IndividualLevel = allCanidsTrueROH %>%
  filter(ROH_length >= 10^6) %>%
  group_by(INDV) %>% 
  summarise(TotalROH = sum(ROH_length)) %>%
  mutate(TotalROHMb = TotalROH/10^6, 
         FROH = TotalROH/2500000000, 
         Population = substr(INDV,1,2))

IndividualsSampled = allCanidsTrueROH %>%
  distinct(INDV) %>% 
  mutate(Population = substr(INDV, 1, 2)) %>%
  group_by(Population) %>%
  tally()

PopulationLevel = allCanidsTrueROH %>% 
  filter(ROH_length >= 10^6) %>%
  group_by(Population) %>% 
  summarise(TotalROH = sum(as.numeric(ROH_length))) %>%
  mutate(TotalROHMb = TotalROH/10^6, 
         FROH = TotalROH/(2500000000*IndividualsSampled$n[match(Population, IndividualsSampled$Population)]))

#FROH unfiltered
unfiltered_PopulationLevel = nonQCROH_allCanids %>%
  filter(ROH_length >= 10^6) %>%
  group_by(Population) %>% 
  summarise(TotalROH = sum(as.numeric(ROH_length))) %>%
  mutate(TotalROHMb = TotalROH/10^6, 
         FROH = TotalROH/(2500000000*IndividualsSampled$n[match(Population, IndividualsSampled$Population)]))

#Plot FROH with ROH greater than 1Mb
orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf", "Isle Royale")

IndividualLevel$Population = mgsub(IndividualLevel$Population, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))

IndividualLevel$Population = factor(IndividualLevel$Population, levels = orderPops)

cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

FROHplot = ggplot(IndividualLevel, aes(x=Population, y=FROH, colour=Population)) + 
  scale_colour_manual(name = "Population", values = cbPalette) + 
  geom_violin(size=1) + 
  geom_point(size=0.5, show.legend=F) + 
  theme_bw() + 
  labs(x = "Population", y=expression(F[ROH])) + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))


#Generate ROH Ranges
#Split into length intervals, find length for each indiv, add column in Gb, add range column
aggIntROH = allCanidsTrueROH %>% 
  select(ROH_length, INDV) %>% 
  filter(ROH_length >= 100000 & ROH_length < 1000000) %>% 
  group_by(INDV) %>% 
  summarise(totalLen = sum(ROH_length)) %>% 
  mutate(totalLenGb = totalLen/10^9, 
         Range = "[0.1-1)Mb") %>% 
  as.data.frame()

aggLongROH = allCanidsTrueROH %>% 
  select(ROH_length, INDV) %>% 
  filter(ROH_length >= 1000000 & ROH_length < 10000000) %>% 
  group_by(INDV) %>% 
  summarise(totalLen = sum(ROH_length)) %>% 
  mutate(totalLenGb = totalLen/10^9, 
         Range = "[1-10)Mb")  %>% 
  as.data.frame()

aggVLongROH = allCanidsTrueROH %>% 
  select(ROH_length, INDV) %>% 
  filter(ROH_length >= 10000000) %>% 
  group_by(INDV) %>% 
  summarise(totalLen = sum(ROH_length)) %>% 
  mutate(totalLenGb = totalLen/10^9, 
         Range = "[10-63)Mb") %>% 
  as.data.frame()

#Merge interval data frames
FinalDF = rbind(aggIntROH, aggLongROH, aggVLongROH)
FinalDF$Range = factor(FinalDF$Range, levels=c("[10-63)Mb","[1-10)Mb","[0.1-1)Mb"))
FinalDF$Population = substr(FinalDF$INDV, 1,2)
FinalDF$Population = mgsub(FinalDF$Population, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))
FinalDF$Population = factor(FinalDF$Population, levels = orderPops)
FinalDF$INDV = factor(FinalDF$INDV, levels = Individuals$V1)

#Function to plot ROH ranges with all ROH greater than 100Kb
plotROHRanges = function(dataFrame){
  ggplot(dataFrame, aes(x=INDV, y=totalLen/10^9, fill=Range)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  coord_flip() + 
  scale_fill_manual(values = c("[0.1-1)Mb"= "bisque3", 
                               "[1-10)Mb" = "darkgoldenrod",
                               "[10-63)Mb" = "indianred4"), 
                    breaks = c("[0.1-1)Mb",
                               "[1-10)Mb", 
                               "[10-63)Mb"), 
                    name = "Range") + 
  labs(x ="Individual", y = "ROH length per bin (Gb)") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10),
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))
  }


#ROH Ranges plots
ROHRangePlot = plotROHRanges(FinalDF) #plot all individuals
ROHRangePlot_TM = plotROHRanges(FinalDF %>% filter(Population == "Tibetan Mastiff")) + ylim(0,1.3)
ROHRangePlot_LB = plotROHRanges(FinalDF %>% filter(Population == "Labrador Retriever")) + ylim(0,1.3)
ROHRangePlot_BC = plotROHRanges(FinalDF %>% filter(Population == "Border Collie")) + ylim(0,1.3)
ROHRangePlot_PG = plotROHRanges(FinalDF %>% filter(Population == "Pug")) + ylim(0,1.3)
ROHRangePlot_AW = plotROHRanges(FinalDF %>% filter(Population == "Arctic Wolf")) + ylim(0,1.3)
ROHRangePlot_IR = plotROHRanges(FinalDF %>% filter(Population == "Isle Royale")) + ylim(0,1.3)
ROHRangePlot_EW = plotROHRanges(FinalDF %>% filter(Population == "Ethiopian Wolf")) + ylim(0,1.3)
