#Load libraries
library(ggplot2)
library(ggtree)
library(phangorn)
library(TreeTools)
library(ggpubr)

#####Read the nexus format PRDM9 mafft aligned data in ######
newickDat = ReadAsPhyDat("~/Documents/DogProject_Clare/LocalRscripts/PRDM9/NexusFormat_PRDM9_rename.txt")
phyDatTree = phyDat(newickDat, type="DNA", levels = NULL)

#####Check out models for distance matrix based on loglikelihood
#model = modelTest(phyDatTree, multicore=TRUE, mc.cores=2)
#model[which.max(model$logLik),] 
#####"env" is an environment which contains all the trees, the data and the calls to allow get the estimated models, e.g. as a starting point for further analysis 
#env = attr(model, "env")
#ls(env=env)
#(F81 = get("F81", env)) #check out F81
#eval(F81, env=env)

#Make neighbor joining tree
dna_dist = dist.ml(phyDatTree, model="JC69")
tree_NJ  = NJ(dna_dist)

#Maximum likelihood 
fit = pml(tree_NJ, phyDatTree)
#fit = pml(tree = rtree(n = length(phyDatTree), tip.label = names(phyDatTree), rooted = FALSE), data = phyDatTree) #fit for random trees 
fit.jc = optim.pml(fit, optNni = TRUE, optEdge = TRUE, model = "JC", rearrangement = "stochastic")#optimize topology and branch lengths

#Bootstrap the tree
bootstrap.jc = bootstrap.pml(fit.jc, bs=100, optNni=TRUE, optEdge = TRUE, multicore=TRUE, control = pml.control(trace=0))

#Compare NJ tree to a tree based on parsimony
#parsimony_optim <- pratchet(phyDatTree)
#plot.phylo(parsimony_optim, use.edge.length=TRUE, cex=0.75)

#Plot all bootstraps
#ggtree(bootstrap.jc) + 
# facet_wrap(~.id, ncol=10) 

#Plot high confidence bootstraps with ggtree
bootstrap.jc.highconf= midpoint(fit.jc$tree)

PRDM9Tree = ggtree(bootstrap.jc.highconf, layout = "slanted") + 
  theme_tree2() + 
  geom_tiplab(size=14, color="blue") + 
  ggtitle("PRDM9") + 
  theme(axis.text.x = element_text(size= 20), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=24, hjust=0.5), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

#######Read the GAPDH nexus format mafft aligned data in ########
newickDat = ReadAsPhyDat("~/Documents/DogProject_Clare/LocalRscripts/PRDM9/NexusFormat_GAPDH_rename.txt")
phyDatTree = phyDat(newickDat, type="DNA", levels = NULL)

#Make neighbor joining tree
dna_dist = dist.ml(phyDatTree, model="JC69")
tree_NJ  = NJ(dna_dist)

#Maximum likelihood 
fit = pml(tree_NJ, phyDatTree)
fit.jc = optim.pml(fit, optNni = TRUE, optEdge = TRUE, model = "JC", rearrangement = "stochastic")#optimize topology and branch lengths

#Bootstrap the tree
bootstrap.jc = bootstrap.pml(fit.jc, bs=100, optNni=TRUE, optEdge = TRUE, multicore=TRUE, control = pml.control(trace=0))

#Plot high confidence bootstraps with ggtree
bootstrap.jc.highconf= midpoint(fit.jc$tree)

GAPDHTree = ggtree(bootstrap.jc.highconf, layout = "slanted") + 
  theme_tree2() + 
  geom_tiplab(size=14, color="blue") + 
  ggtitle("GAPDH") + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=24, hjust=0.5), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

geneTrees = ggarrange(PRDM9Tree, 
                      GAPDHTree, 
                      align = 'hv', 
                      nrow = 2,
                      labels=c("A","B"))
