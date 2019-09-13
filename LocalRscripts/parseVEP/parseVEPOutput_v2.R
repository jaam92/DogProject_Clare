#Load Libraries
library(dplyr)
library(data.table)
library(hablar)
library(tidyr)
library(mgsub)

#Load files
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/VEP/parseVEP")

#zero and four fold files are from snpEFF annots from Eduardo
zeroFold = fread("0fold.deg", header = F) %>% 
  select(V1) %>%
  mutate(annot = "zeroFold") %>%
  as.data.frame()
head(zeroFold)

fourFold = fread("4fold.deg", header = F) %>%
  select(V1) %>%
  mutate(annot = "fourFold") %>%
  as.data.frame() 
head(fourFold)

snpEffAnnot = rbind.data.frame(zeroFold,fourFold) 
names(snpEffAnnot)[1] = "chrom_pos"
rm(zeroFold)#delete
rm(fourFold)#delete

#load filenames for VEP and SIFT Annots
filenames = list.files(pattern = glob2rx("masterFile_allCanidsN75_ANmt135_Chr*_fixHeader.txt"))
print(filenames)

#dataFrame with all annots
allChroms = data.frame()

#Loop thru all VEP files
for (i in 1:length(filenames)) {
  
df = read.delim(filenames[i]) %>%
    retype() %>%
    filter(IMPACT != "MODIFIER") %>% #remove intergenic stuff
    mutate(IMPACT = mgsub(IMPACT, c("LOW", "MODERATE", "HIGH"), c("0", "1", "2"), perl = TRUE), #make heirarchy numeric so we can select the most damaging annotation for a given transcript
           CHROM_POS = gsub(":", "_", X.Location, perl = TRUE)) %>% #change sep so I can compare to Eduardo's annotation
    select(-c(X.Location))#remove location col since we no longer need it
  
  #Keep the most damaging annotation
  keepHighImpact =  df %>%
    group_by(CHROM_POS) %>% 
    slice(which.max(as.numeric(IMPACT))) %>% #choose most damagind annot and if there is a tie for impact pick first reported (normally ends up being the canonical)
    mutate(IMPACT = mgsub(IMPACT, c("0", "1", "2"), c("LOW", "MODERATE", "HIGH"),perl = TRUE), #change the impact column back to original annot
       SNPEFF = snpEffAnnot$annot[match(CHROM_POS,snpEffAnnot$chrom_pos)], #grab annotation from Eduardo's snpEFF file
       SIFT_SCORE = suppressWarnings(as.numeric(gsub("-", "NA", SIFT, perl = TRUE))), #make sift score numeric NAs get introduced so I supress the warning
       SIFT_CONSEQUENCE = ifelse(as.numeric(SIFT) > 0.5, "benign", "deleterious")) %>% #classify missense mutations as benign or damaging based on sift score
    separate(CHROM_POS, c("CHROM", "POS"), sep = "_") %>% 
    select(CHROM, POS, Consequence, IMPACT, SIFT_SCORE, SIFT_CONSEQUENCE, SYMBOL, Gene, CANONICAL, SNPEFF) %>% #keep columns I want
    plyr::rename(c("Consequence"="CONSEQUENCE", "Gene"="GENE")) #rename columns
  
  
  #Generate a table with all chroms
  allChroms = rbind(allChroms, keepHighImpact)#Generate a table that contains all annotation data
  
  #Generate output file per chrom
  #write.table(keepHighImpact, file=paste("FinalAnnot_VEPSIFTSnpEff_", filenames[i], sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

print("finished")
  
}

sorted_allChroms = allChroms %>% 
  mutate(CHROM_num = as.numeric(gsub("chr","",CHROM))) %>%
  arrange(CHROM_num) %>% #order by CHROM
  select(-c(CHROM_num))

#Generate output file all chroms
write.table(sorted_allChroms, file= "FinalAnnot_VEPSIFTSnpEff_allChroms.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#rm(sorted_allChroms) #delete we no longer need this

#compare annotations
#zeroFoldAnnot = allChroms %>% filter(SNPEFF == "zeroFold")


#Compare number of annotations VEP and snpEFF
#VEP = allChroms %>% 
#  group_by(CHROM) %>% 
#  count() %>%
#  mutate(method = "VEP")

#snpEff = snpEffAnnot %>% 
#  mutate(CHROM = gsub("_.*","", chrom_pos)) %>% 
#  group_by(CHROM) %>% 
#  count()%>%
#  mutate(method = "snpEff")

#ggplot() + geom_bar(data=VEP, aes(x=CHROM,y=n), stat = "identity", fill="red") + geom_bar(data=snpEff, aes(x=CHROM,y=n), stat = "identity", fill="blue", alpha = 0.2) + theme_bw()
  
  

