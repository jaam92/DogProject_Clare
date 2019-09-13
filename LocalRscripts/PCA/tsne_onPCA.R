#pcair
canidMatrix = pc.df %>%
  select(-c(sample.id, Population)) %>%
  as.matrix()

#smartPCA
#canidMatrix = pc.df %>%
#  select(-c(Individual, Population)) %>%
#  as.matrix()

#run tsne
###I am supplying the PCA so I set the pca argument to F and set perplexity to 24 bc (3*perplexity < nrow(df)-1)###
set.seed(404) #seed for pcair
#set.seed(308) #seed for smart pca
tsne_Canids = Rtsne(canidMatrix, 
                    pca=F,
                    normalize = F,
                    perplexity = 24,
                    max_iter = 5000)

#Make data frame
plotTsne = tsne_Canids$Y %>%
  as.data.frame() %>%
  mutate(Individual = pc.df$Individual,
         Population = pc.df$Population)

plotTsne$Population = mgsub(plotTsne$Population, pattern =c("BC", "LB", "PG", "TM", "AW", "EW", "IR"), replacement =c("Border Collie", "Labrador Retriever", "Pug", "Tibetan Mastiff", "Arctic Wolf",  "Ethiopian Wolf", "Isle Royale"))

#plot tsne
cbPalette = c("Arctic Wolf" = "gray25", "Ethiopian Wolf" = "#D55E00",  "Isle Royale" = "steelblue", "Border Collie" = "#009E73", "Labrador Retriever" = "gold3", "Pug" = "mediumpurple4", "Tibetan Mastiff" = "#CC79A7")

tsne1VStsne2 = ggplot(plotTsne, aes(x=V1, y=V2, color=Population)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Population") + 
  theme_bw() + 
  labs(x="tsne1",y="tsne2") + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
