#Load Libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(mgsub)
library(data.table)

#Functions
###Make a folded SFS###
#THIS IS BERNARD's ORIGINAL >> need to add 0 to front and end to stand in for monomorphic (don't need to do this if it came from dadi)
fold <- function(SFSCountCol, n, norm=TRUE){
  if (length(SFSCountCol) < (n+1)){
    data = c(SFSCountCol, rep(0,(n+1)-length(SFSCountCol)))
  }
  data = SFSCountCol[2:n] # adds together sfs and backwards sfs
  data_fold = data + rev(data) # takes first half of that added together sfs (but not middle entry)
  data_fold = data_fold[1:(n/2-1)]# adds middle entry that didn't have anything added to the end
  data_fold = c(data_fold,data[(n/2)])# truncates and sums up remaining fields if desired (not needed here)
  #data_trunc = c(data_fold[1:(trunc-1)],sum(data_fold[trunc:length(data_fold)]))
  #if (norm){
  #  data_trunc = data_trunc/sum(data_trunc)
  #}
  #return(data_trunc)
  return(data_fold)
}
#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/SFS")
fnames = list.files(pattern = "\\_SFS_allChroms_SummaryFile_N6.txt$")
#fnames = list.files(pattern = "\\_NeutralSFS_allChroms_SummaryFile_N6.txt$")
df = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FileName")#Read all the files and create a FileName column to store filenames

#Plot Data
PlottingUnfolded = df %>%
  filter(FreqBin > 0 & FreqBin < 12) %>% #remove fixed stuff
  mutate(Population = substr(FileName,1,2)) %>%
  group_by(Population, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  mutate(Population = mgsub(Population,pattern=c("BC", "LB", "PG", "TM", "AW", "EW", "IR"),replacement = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf", "Isle Royale")))

TotalCounts = PlottingUnfolded %>%
  group_by(Population) %>%
  summarise(TotalSites = sum(WholeGenomeCounts))

PlottingUnfolded$Proportional = PlottingUnfolded$WholeGenomeCounts/(TotalCounts$TotalSites[match(PlottingUnfolded$Population, TotalCounts$Population)])

#Prep for plotting
orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf", "Isle Royale")

PlottingUnfolded$Population = factor(PlottingUnfolded$Population, levels = orderPops)

cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

#Plot only the bins starting with singletons
UnfoldedSFS = ggplot(PlottingUnfolded, aes(y=Proportional, x=FreqBin, fill=Population)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:11) + 
  scale_fill_manual(values= cbPalette) + 
  labs(x= "SNP Frequency", y= "Proportion of SNPs", title = "All Populations Whole Genome SFS") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        plot.title=element_text(size=32, face = "bold", hjust=0.5), 
        axis.title=element_text(size=32), 
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24))

#print(UnfoldedSFS)

Unfolded = df %>%
  mutate(Population = substr(FileName,1,2)) %>%
  group_by(Population, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  mutate(Population = mgsub(Population,pattern=c("BC", "LB", "PG", "TM", "AW", "EW", "IR"),replacement = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf", "Ethiopian Wolf", "Isle Royale"))) 

splitUnfolded = split.data.frame(Unfolded, Unfolded$Population)

PlottingFolded = sapply(names(splitUnfolded), function(x) fold(splitUnfolded[[x]]$WholeGenomeCounts, 12)) %>%
  as.data.frame() %>%
  mutate(bins = 1:6) %>%
  melt(id=c("bins")) 

TotalCountsFolded = PlottingFolded %>%
  group_by(variable) %>%
  summarise(TotalSites = sum(value))

PlottingFolded$Proportional = PlottingFolded$value/(TotalCountsFolded$TotalSites[match(PlottingFolded$variable, TotalCountsFolded$variable)])

FoldedSFS = ggplot(PlottingFolded, aes(y=Proportional, x=bins, fill=variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_fill_manual(values= cbPalette, name = "Population") + 
  scale_x_continuous(breaks=1:6) +
  labs(x= "SNP Frequency", y= "Proportion of SNPs", title= "All Populations Whole Genome SFS") + 
  theme_bw()  + 
  theme(axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        plot.title=element_text(size=32, face = "bold", hjust=0.5), 
        axis.title=element_text(size=32), 
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24))

#print(FoldedSFS)

#SFS with only Wolves
WolfOnlySFS = ggplot(PlottingFolded %>% 
         filter(variable == "Arctic Wolf"| variable == "Ethiopian Wolf"| variable == "Isle Royale" ), aes(y=Proportional, x=bins, fill=variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:6) +
  scale_fill_manual(values= cbPalette, name = "Population") + 
  labs(x= "SNP Frequency", y= "Proportion of SNPs", title="Wolves Whole Genome SFS") + 
  theme_bw()  + 
  theme(axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        plot.title=element_text(size=32, face = "bold", hjust=0.5), 
        axis.title=element_text(size=32), 
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24))

#print(WolfOnlySFS)