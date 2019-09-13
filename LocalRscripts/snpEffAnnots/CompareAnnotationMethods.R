library(dplyr)
library(data.table)
library(ggplot2)
library(randomcoloR)

#Compare my snpEff Annots and my VEP annots
setwd("~/Documents/DogProject_Clare/LocalRscripts/snpEffAnnots")
snpEff = read.delim("~/Documents/DogProject_Clare/LocalRscripts/snpEffAnnots/anns_snpEff_cannonicalOnly_reformat.txt")
VEP = read.delim("~/Documents/DogProject_Clare/LocalRscripts/parseVEP/FinalAnnot_VEPSIFTSnpEff_allChroms.txt")

#Differences 
snpEff$CHROM_POS = paste(snpEff$chr, "_", snpEff$pos, sep = "")
VEP$CHROM_POS = paste(VEP$CHROM, "_", VEP$POS, sep = "")
difference = dim(VEP)[1] - dim(snpEff)[1]
print(paste("There are this many more annotations in VEP:",difference, sep=" "))

#Shared annotations
exclusivelyVEP = anti_join(VEP, snpEff, by="CHROM_POS")
VEPDiff = dim(exclusivelyVEP)[1]
print(paste("There are this many exclusive annotations in VEP:",VEPDiff, sep=" "))

exclusivelysnpEff = anti_join(snpEff, VEP, by="CHROM_POS")
snpEffDiff = dim(exclusivelysnpEff)[1]
print(paste("There are this many exclusive annotations in snpEFF:",snpEffDiff, sep=" "))

#Plot extra annotations from VEP
plotDifferences = exclusivelyVEP %>%
  group_by(CONSEQUENCE) %>%
  count()

legendPal = distinctColorPalette(dim(plotDifferences)[1])
  
ggplot(plotDifferences, aes(x=reorder(CONSEQUENCE, -n),y=n, fill=CONSEQUENCE)) + geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + scale_fill_manual(values = legendPal) + theme_bw() + labs(x="Annotation", y="Count") + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=18), legend.text=element_text(size=14), legend.position = "bottom") + guides(fill=guide_legend(nrow=8,byrow=TRUE))

