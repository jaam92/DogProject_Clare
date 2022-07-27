library(tidyverse)
library(mgsub)

#Set working directory
setwd("~/Documents/DogProject_Clare/LocalRscripts/pairwiseFst")

#Set breeds as factor so wolves and dogs group together
orderPops = c("border collie", "labrador retriever", "pug", "tibetan mastiff", "Arctic wolf", "Isle Royale wolf", "Ethiopian wolf")

#read file in
df = read_delim("~/Documents/DogProject_Clare/LocalRscripts/pairwiseFst/pairwiseFst_reformat.out", delim = "\t") %>%
  mutate(B1 = mgsub(B1, pattern=c("BC", "LB", "PG", "TM", "AW", "EW", "IR"),replacement = c("border collie", "labrador retriever", "pug", "tibetan mastiff", "Arctic wolf", "Ethiopian wolf", "Isle Royale wolf")),
         B2 = mgsub(B2, pattern=c("BC", "LB", "PG", "TM", "AW", "EW", "IR"),replacement = c("border collie", "labrador retriever", "pug", "tibetan mastiff", "Arctic wolf", "Ethiopian wolf", "Isle Royale wolf")),
         B1 = factor(B1, level=orderPops),
         B2 = factor(B2, level=orderPops)) %>%
  mutate_if(is.numeric, round, digits=3) %>% 
  na.omit() #%>%
  #mutate(WeightedFst = ifelse(WeightedFst==0, NA, WeightedFst))

#plot
ggplot(df, aes(x = B1, y = B2, fill = WeightedFst)) +
  geom_tile() +
  scale_fill_distiller(palette = "GnBu", name= expression('F'[ST]), limits=c(0,1), na.value = "white") +
  geom_text(aes(label = ifelse(WeightedFst != 1, WeightedFst, "")), size = 10) + 
  theme_bw() + 
  labs(x= NULL, y = NULL) + 
  theme(axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        plot.title=element_text(size=32, face = "bold", hjust=0.5), 
        axis.title=element_text(size=32), 
        legend.title=element_text(size=32), 
        legend.text=element_text(size=30),
        panel.grid = element_blank())
