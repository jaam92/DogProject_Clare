library(ape)
library(igraph)
library(tidyverse)
library(mgsub)
library(corrplot)

#read file in
setwd("~/Documents/DogProject_Clare/LocalRscripts/pairwiseFst")
df = read.table("~/Documents/DogProject_Clare/LocalRscripts/pairwiseFst/pairwiseFst_reformat.out", header = T, stringsAsFactors=FALSE) %>%
  mutate(B1 = mgsub(B1, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale")),
         B2 = mgsub(B2, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale")))

#make data framewith comparisons the other way from current
df2 = df %>% 
  select(B1 = B2, B2 = B1, WeightedFst)

#merge together
mergeDF = rbind.data.frame(df, df2, stringsAsFactors=FALSE)  

#Set breeds as factor so Wolves and dogs group together
orderPops = c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale")

mergeDF$B1 = factor(mergeDF$B1, levels=orderPops)
mergeDF$B2 = factor(mergeDF$B2, levels=orderPops)

#plot
#ggplot(data = mergeDF, aes(x=B1, y=B2,fill=WeightedFst), colour = "white") +  geom_tile() + scale_fill_gradient(low = "blue", high = "red") + theme_bw() + theme(axis.text.x = element_text(size  = 24,angle=40, vjust=1, hjust=1), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18)) + labs(x="Group 1", y="Group 2")


#Plot lower triangle as heat map
g = graph.data.frame(df, directed=FALSE)
fstMatrix=get.adjacency(g, attr="WeightedFst", sparse=FALSE)
#write.csv(fstMatrix,file="FstMatrix.csv")
col4 = colorRampPalette(c("cyan",  "yellow", "red"))

corrplot(fstMatrix, type="lower", order="hclust", tl.col="black", tl.srt=45, cl.lim = c(0,1), method = "number", col = col4(20))
