##This script will grab Exons from ensembl
###Generated Ensembl data by downloading exons from UCSC as bed file using protocol from Tanya(https://github.com/tanyaphung/SexBiased/tree/master/filter_nucleotide_sites)
#Load Libraries
library(dplyr)
library(tidyr)

####Read files in
setwd("~/Documents/DogProject_Clare/LocalRscripts/generateExonCoords")
genes = read.delim("~/Documents/DogProject_Clare/LocalRscripts/generateExonCoords/EnsemblGenes_CanFam3.1.bed", check.names = F, stringsAsFactors = F) #come from UCSC table browser (downloaded as bed)
gene_names = read.table("~/Documents/DogProject_Clare/LocalRscripts/generateExonCoords/EnsemblGenes_CanFam3.1_geneNames.txt", stringsAsFactors = F)

###Add gene name -> calculate transcript length -> keep only longest transcript -> keep autosomal -> parse exons
GeneSet = genes %>% 
  mutate(chrom = gsub("chr", "", chrom), 
         AbbrevName = gene_names$V2[match(name, gene_names$V1)],
         transcript_length = transcriptionEnd - transcriptionStart) %>% 
  group_by(AbbrevName) %>% 
  filter(transcript_length == max(transcript_length) & chrom %in% 1:38) %>% 
  ungroup() %>%
  distinct(AbbrevName,.keep_all= TRUE) %>% #removes duplicates of transcripts that are the same length
  separate_rows(exonStarts, exonEnds) %>% 
  mutate(CHROM = sub("^", "chr", chrom)) %>%
  select(CHROM, transcriptionStart, transcriptionEnd, transcript_length, codingRegionStart, codingRegionEnd, exonStarts, exonEnds, name, AbbrevName) %>%
  as.data.frame()

FinalGeneSet = GeneSet[!(GeneSet$exonStarts == ""),] #remove empty rows generated from separate rows

write.table(FinalGeneSet, "EnsemblGenes_CanFam3.1_SingleTranscript.bed", sep = "\t", quote = F, row.names = F, col.names = T)
