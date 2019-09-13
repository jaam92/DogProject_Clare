####This RScript will compare my VEP annotations with Eduardo's SNPEff annots that had at most 3 transcripts

 
library(dplyr)
library(data.table)

#Load Files 
EduardoDF = fread("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/VEP/parseVEP/SharedAnnotsEduardo/3annotsEduardo.txt")
myDF = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/VEP/parseVEP/FinalAnnot_VEPSIFTSnpEff_allChroms.txt")

#Add col with chrom and pos to compare
myDF$CHROM_POS = paste(myDF$CHROM, "_" , myDF$POS, sep="")

#Compare
dim(EduardoDF)
dim(myDF) #76245 annots
dim(EduardoDF %>% filter(V1 %in% myDF$CHROM_POS)) #47032 shared annots

#Check out shared annots
sharedAnnotsEduardo = EduardoDF %>% 
filter(V1 %in% myDF$CHROM_POS) %>% 
group_by(V3,V4,V5) %>% 
count()

colnames(sharedAnnotsEduardo) = c("ANNOT1","ANNOT2", "ANNOT3", "COUNT")

#Output with my annotation
addMyAnnotsEduardo = EduardoDF %>%
filter(V1 %in% myDF$CHROM_POS) %>%
mutate(myANNOT = myDF$CONSEQUENCE[match(V1, myDF$CHROM_POS)]) %>%
select(V1, V3, V4, V5, myANNOT)
colnames(addMyAnnotsEduardo) = c("CHROM_POS","ANNOT1","ANNOT2", "ANNOT3", "myANNOT")

#write to output file
write.table(addMyAnnotsEduardo, file="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/VEP/parseVEP/SharedAnnotsEduardo/BreakDownSharedAnnot.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(sharedAnnotsEduardo, file="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/VEP/parseVEP/SharedAnnotsEduardo/CountsBreakDownSharedAnnot.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
