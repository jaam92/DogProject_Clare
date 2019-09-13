#Load libraries
library(data.table)
library(ggplot2)
library(hablar)
library(dplyr)
library(mgsub)

#Load files
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/vcfToolsROH/AnnotateROH/AnnotatedRange")

#Load vcftools outputs
fnames = list.files(pattern = "\\ROHQC.LROH$")
df = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE, idcol = "FileName")#Read all the files and create a FileName column to store filenames
nonQCROH_allCanids = subset(df, df$CHROM != "CHROM") %>% retype() %>% mutate(Population = substr(INDV, 1, 2)) %>% as.data.frame()#remove header lines and convert character columns to numeric
rm(df)#delete old df

#QC ROH Segments
#Old Filter 
##Will remove ROH if the proportion of ROH segment covered by snps is not within a SD of the mean
#z = data.table(nonQCROH_allCanids)
#z[,ToKeep := abs(nonQCROH_allCanids$Prop_Cov - mean(nonQCROH_allCanids$Prop_Cov)) < sd(nonQCROH_allCanids$Prop_Cov)][ToKeep  == TRUE] #create variable that identifies what to drop or keep 
#table(z$ToKeep, z$Population)/length(z$ToKeep) #Proportion that fall within drop or keep for each Population
#allCanidsTrueROH = subset(z, z$ToKeep == "TRUE")

#Filter chosen for project (ROH is at least 10kb)
allCanidsTrueROH = subset(nonQCROH_allCanids, nonQCROH_allCanids$Prop_Covered >= 0.4 & nonQCROH_allCanids$ROH_length != 1 & nonQCROH_allCanids$ROH_length != 0 & nonQCROH_allCanids$ROH_length >= 10000) %>% 
 # select(-c(FileName)) %>%
  #mutate(CHROM = paste("chr",CHROM,sep = ""))

write.table(allCanidsTrueROH, file="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/vcfToolsROH/AnnotateROH/AnnotatedRange/FinalROHQCd_min10Kb_allIndivs_Dec2018_DogProjClare.LROH", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


#Filter to ROH of at least 1Mb or greater
allCanidsTrueROH = subset(nonQCROH_allCanids, nonQCROH_allCanids$Prop_Covered >= 0.4 & nonQCROH_allCanids$ROH_length != 1 & nonQCROH_allCanids$ROH_length != 0 & nonQCROH_allCanids$ROH_length >= 10^6) %>%
  select(-c(FileName)) %>%
  mutate(CHROM = paste("chr",CHROM,sep = ""))

write.table(allCanidsTrueROH, file="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/vcfToolsROH/AnnotateROH/AnnotatedRange/FinalROHQCd_min1Mb_allIndivs_Dec2018_DogProjClare.LROH", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

