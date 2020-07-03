#Load necessary packages
library(tidyverse)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(corrplot)
library(ggdendro)
library(dendextend)
library(seriation) #olo function
library(cluster) #use silhouette

#Set working directory and load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/PCA")

####SNPRelate to run KING####
gds = snpgdsOpen("masterFile_allChroms_updateFID.gds")
sampIds = read.gdsn(index.gdsn(gds, "sample.id")) #grab sample ids 
famIds = substr(sampIds,1,2) #make family ids

#Run KING
KING = snpgdsIBDKING(gds, sample.id = sampIds, autosome.only = F, num.thread = 2)

#IBS matrix
#ibs.hc = snpgdsHCluster(snpgdsIBS(gds, num.thread=2))

snpgdsClose(gds)

#Make the IBS0 (probability IBS=0) matrix
kingIBSMatrix = KING$IBS0
colnames(kingIBSMatrix) <- rownames(kingIBSMatrix) <- KING$sample.id
disMatrix = as.dist(kingIBSMatrix) #convert to distance matrix

#Use all methods to cluster
#code from https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html
#hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
#king_dendlist = dendlist()
#for(i in seq_along(hclust_methods)) {
#  hc_king = hclust(disMatrix, method = hclust_methods[i])   
#  king_dendlist = dendlist(king_dendlist, as.dendrogram(hc_king))
#}
#names(king_dendlist) = hclust_methods
#king_dendlist_cor = cor.dendlist(king_dendlist)#Look at similarity
#corrplot(king_dendlist_cor, "pie", "lower", title = "compare clustering methods") #All methods produce similar tree 

#Perform heirarchical clustering
hc_avg = hclust(disMatrix, method = "average")
hc_avg_reorder = reorder(hc_avg, disMatrix, method = "OLO") #reorder with optimal leaf ordering

#Check whether cluster membership is well supported (Sw ~ 0.6, so it's pretty good)
####For interpretation of silhouette width https://www.stat.berkeley.edu/~spector/s133/Clus.html
####Code: https://www.rdocumentation.org/packages/cluster/versions/2.0.8/topics/silhouette

#c6 = c("tomato", "darkseagreen4", "steelblue", "purple2", "goldenrod4", "gray25")
#par(mfrow= c(3,2), oma= c(0,0, 3, 0),mgp= c(1.6,.8,0), mar= .1+c(4,2,2,2))
#for(k in 2:6){
#plot(silhouette(cutree(hc_avg_reorder, k= k),disMatrix), 
#     main = paste("k = ",k), 
#     do.n.k=FALSE, 
#     col = c6[1:k])
#}


#Function to plot 
dend = hc_avg_reorder %>% 
    as.dendrogram %>%
    set("branches_k_color", k = 4, c("#D55E00","steelblue","gray25","darkseagreen4")) %>% 
    set("branches_lwd", 0.7) %>%
    set("labels_cex", 1.5) %>% 
    set("labels_colors", k = 4, c("#D55E00","steelblue","gray25","darkseagreen4")) %>%
    set("leaves_pch", 19) %>% 
    set("leaves_cex", 0.5) 

ggd1 = as.ggdend(dend)
  
IBSTree = ggplot(ggd1, horiz = F) + 
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 24),
          axis.title=element_text(size=24))


#pvclust to give clusters with high confidence
###info:http://stat.sys.i.kyoto-u.ac.jp/prog/pvclust/
#library(pvclust) 
#x = pvclust(kingIBSMatrix, method.dist = "cor", method.hclust = "average", nboot = 10000)
#plot(x)
#pvrect(x, alpha = 0.95)


#IBS distance matrix from snpgds 
#rv = snpgdsCutTree(ibs.hc)
#IBSTreeSNPgds = ggplot(rv$dendrogram, horiz = F) + 
#ggtitle("Dist Matrix") + 
#  theme(axis.text.x = element_blank(), 
#        axis.text.y = element_text(size  = 20),
#        axis.title=element_text(size=24))