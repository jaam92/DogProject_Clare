###Warning just appear because I am cutting out the column information I do not need when I use separate rows

#Load Libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(tidyr)

#Load files
setwd("~/Documents/DogProject_Clare/LocalRscripts/parseVEP/")

#zero and four fold files are from snpEFF annots from Eduardo
zeroFold = fread("0fold.deg", header = F) %>% 
  select(V1) %>%
  mutate(annot = "zeroFold") %>%
  as.data.frame()

fourFold = fread("4fold.deg", header = F) %>%
  select(V1) %>%
  mutate(annot = "fourFold") %>%
  as.data.frame() 

snpEffAnnot = rbind.data.frame(zeroFold,fourFold) 
names(snpEffAnnot)[1] = "chrom_pos"
rm(zeroFold)#delete
rm(fourFold)#delete

#load filenames for VEP and SIFT Annots
filenames = list.files(pattern = glob2rx("chr*_codingAnnots.txt"))

#dataFrame with all annots
allChroms = data.frame()

#Loop thru all VEP files
for (i in 1:length(filenames)) {
      
      #read file
      df = read.table(filenames[i])
      
      reformattedDF = df %>% 
        mutate(CHROM_POS = paste(V1, V2, sep="_"), #add col to annotate with snpEFF
                    V3 = gsub("^[^|]*|", "", V3), #remove all the stuff before VEP annot
                    V3 = gsub("\\|", " ", V3), #replace all the | with a space
                    CANONICAL = str_extract(V3, pattern="YES"), #check whether canonical annot
                    SIFT = gsub(".*YES", "", V3), #grab SIFT score
                    SIFT = gsub(",.*", "", SIFT), #remove all the stuff after the score
                    V3 = gsub("Trans.*","", V3), #keep info for the first (most damaging) annot
                    SNPEFF = snpEffAnnot$annot[match(CHROM_POS,snpEffAnnot$chrom_pos)]) %>% 
        separate(V3, c("trash", "ANNOT", "IMPACT", "GENE", "ENS_GENEID"), sep = " ") %>% 
        select(V1, V2, ANNOT, IMPACT, SIFT, GENE, ENS_GENEID, CANONICAL, SNPEFF) %>% #remove extra cols
        plyr::rename(c("V1"="CHROM", "V2"="POS")) #rename cols
      
      #Table with all chroms
      allChroms = rbind(allChroms, reformattedDF)#Generate a table that contains everyones data
      
      #Generate file for Tanya's script 
      #write.table(reformattedDF, file=paste("parsedVEP_", filenames[i], sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

}

#compare annotations
zeroFoldAnnot = allChroms %>% filter(SNPEFF == "zeroFold")


#Compare number of annotations VEP and snpEFF
VEP = allChroms %>% 
  mutate(CHROM = as.numeric(as.character(gsub("chr", "", CHROM)))) %>%
  group_by(CHROM) %>% 
  count() %>%
  mutate(method = "VEP")

snpEff = snpEffAnnot %>% 
  mutate(CHROM = gsub("_.*","", chrom_pos),
         CHROM = as.numeric(as.character(gsub("chr", "", CHROM)))) %>% 
  group_by(CHROM) %>% 
  count()%>%
  mutate(method = "snpEff")

(compCounts = ggplot() + geom_bar(data=VEP, aes(x=CHROM,y=n), stat = "identity", fill="red") + geom_text(data=VEP, aes(x=CHROM, y=n, label=n), vjust=-1) + geom_bar(data=snpEff, aes(x=CHROM,y=n), stat = "identity", fill="blue", alpha = 0.2) + labs(x="Chromosome", y="Number of Annotations") + theme_bw() + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)))

#Check out the four fold sites
fourFoldComps = allChroms %>%
  filter(SNPEFF == "fourFold")

fourFoldComps %>%
  group_by(IMPACT) %>%
  count()

fourFoldComps %>%
  group_by(ANNOT) %>%
  count()

#Check out the zero fold sites
zeroFoldComps = allChroms %>%
  filter(SNPEFF == "zeroFold")

zeroFoldComps %>%
  group_by(IMPACT) %>%
  count()

zeroFoldComps %>%
  group_by(ANNOT) %>%
  count()
