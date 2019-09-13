#Load Data
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

#Set working directory
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ROHPermutation")

#Load files
Indivs = read.delim("IndivID.txt")

#Set id and group variables for ordering
indNames = Indivs$ID
group = Indivs$Population

#Empty data frames for annotation, snp, and ROH data
annotation_data = data.frame()
snp_data = c()
snp_data_pos = c()
fixedID = list()
ROH_data = c()
ROH_data_pos = c()

####WORK WITH SNP DATA FIRST

#Plot specifications
jpeg(file = "SNP_data.jpg", width=1200*3, height=500*11/2, pointsize = 16)
par(mfrow=c(8,3))
par(mar=c(5,7,2,2)+0.1)

#Grab genotype files sort according to chrom number
fn = dir(pattern = "PermFormatGT_", full.names = TRUE)
fn = str_sort(fn, numeric = TRUE) 

#Loop through all chromosomes 
for (c in fn){
  temp = t(read.table(file = c, header=T, stringsAsFactors = FALSE, check.names = FALSE)) #expecting a file with a header
  
  #Grab info for sorting dataframe
  chromo = str_extract(c,"[[:digit:]]+") #grab chromosome number
  indID = unlist(temp[1,]) #grab IDs from first row
  colnames(temp) = indID #set column names
  
  #Sort data frame according to position and indiv id 
  temp = rbind(temp[1,], (temp[2:nrow(temp),])[order(as.numeric(substr(rownames(temp[2:nrow(temp),]), 1, 100))),])  #Sort according to position if not already sorted
  temp = temp[2:nrow(temp),] #remove the row with the IDs
  
  #Order columns relative to what we see in the indiv file
  if (any(indID != indNames)) {
    print(paste("Order of individuals switched from in file", c))
    temp = temp[,indNames] #sort according to indiv file
    indID = indNames
  }
  
  #Rename GT values in data frame
  temp = plyr::mapvalues(temp, from=c("Missing","AncHom", "Het", "DerHom"), to=c("-1", "0", "1", "2"))
  class(temp) <- "numeric"
  
  #uncomment this section to remove fixed sites (all AncHom or all DerHom) -- need to do this also for the runs later in the script...
  #fixedID[[chromo]] = apply(temp, 1, function(x) all(x==0, na.rm=T) | all(x==2, na.rm=T))
  #temp = temp[!fixedID[[chromo]],]
  
  #grab posititon and set up chromosome variable
  pos = as.numeric(substr(rownames(temp), 1, 100))
  chromosome = rep(chromo, length(pos))
  
  IDorder = order(pos)
  if (any(temp==-1)) {
    colors = c("yellow", "white", "red", "black") 
  } else { 
    colors = c("white", "red", "black")
  }
  
  image(x=pos[IDorder], y=1:ncol(temp), z=as.matrix(temp[IDorder,]), col=colors, xlab="Position", ylab="", xlim=c(0,2.5e8), yaxt='n')
  axis(2, at = 1:ncol(temp), labels=indNames, las = 1)
  title(main=paste("Chromosome:", chromosome[1]))
  rect(xleft=min(pos), xright=max(pos), ybottom=0, ytop=61, col="#00FF0010")
  
  snp_data = rbind(snp_data, temp)
  snp_data_pos = rbind(snp_data_pos, data.frame(chr = chromosome, pos = pos))
}

dev.off() #close file 


####WORK WITH ROH DATA

#plot specifications
jpeg(file = "ROH_data_10kb.jpg", width=1200*3, height=500*11/2, pointsize = 16)
par(mfrow=c(8,3))
par(mar=c(5,7,2,2)+0.1)

#Grab files 
fn = dir(pattern = "PermFormat_TenKbROH_", full.names = TRUE)
fn = str_sort(fn, numeric = TRUE) 

#Loop through all chromosomes 
for (c in fn){
  temp = t(read.table(file = c, header=T, stringsAsFactors = FALSE, check.names = FALSE)) #expecting a file with a header
  
  #Grab info for sorting dataframe
  chromo = str_extract(c,"[[:digit:]]+") #grab chromosome number
  indID = unlist(temp[1,]) #grab IDs from first row
  colnames(temp) = indID #set column names
  
  #Sort data frame according to position and indiv id 
  temp = rbind(temp[1,], (temp[2:nrow(temp),])[order(as.numeric(substr(rownames(temp[2:nrow(temp),]), 1, 100))),])  #Sort according to position if not already sorted
  temp = temp[2:nrow(temp),] #remove the row with the IDs
  
  #Order columns relative to what we see in the indiv file
  if (any(indID != indNames)) {
    print(paste("Order of individuals switched from in file", c))
    temp = temp[,indNames] #sort according to indiv file
    indID = indNames
  }
  
  #Rename GT values in data frame
  temp = plyr::mapvalues(temp, from=c("Missing","nonROH", "ROH"), to=c("-1", "0", "1"))
  class(temp) <- "numeric"
  
  #uncomment to remove fixed sites
  #temp = temp[!fixedID[[chromo]],]
  
  #grab posititon and set up chromosome variable
  pos = as.numeric(substr(rownames(temp), 1, 100))
  chromosome = rep(chromo, length(pos))
  IDorder = order(pos)
  image(x=pos[IDorder], y=1:ncol(temp), z=as.matrix(temp[IDorder,]), col=c("grey", "white", "red"), xlab="Position", ylab="", xlim=c(0,2.5e8), yaxt='n')
  axis(2, at = 1:ncol(temp), labels=indNames, las = 1)
  title(main=paste("Chromosome:", chromosome[1]))
  rect(xleft=min(pos), xright=max(pos), ybottom=0, ytop=61, col="#00FF0010")
  
  ROH_data = rbind(ROH_data, temp)
  ROH_data_pos = rbind(ROH_data_pos, data.frame(chr = chromosome, pos = pos))
}
dev.off()

####WORK WITH ANNOTATIONS

fn = dir(pattern = "_FinalGERPAnnots.txt", full.names = TRUE)
fn = str_sort(fn, numeric = TRUE) 

for (c in fn) {
  temp = read.table(file = c, header=F, stringsAsFactors = FALSE)
  colnames(temp) = c("chrNum", "pos", "annotation")
  
  #Grab info for sorting dataframe
  chromo = str_extract(c,"[[:digit:]]+") #grab chromosome number
  
  #remove fixed sites by comparing to ROH data 
  #warnings pop up but you can ignore
  temp = temp %>% mutate(chr = as.factor(temp$chrNum)) %>%
    right_join(ROH_data_pos %>% 
    filter(chr == chromo), by=c("chr","pos")) %>% 
    select(-c(chr)) %>%
    mutate(chrNum = chromo) %>% 
    arrange(pos) %>% rename("chr" = "chrNum")
  
  #Save new data
  annotation_data = rbind(annotation_data, data.frame(temp))
  
}

# Check that all positions are the same and in same order (both must be FALSE):
any(paste(annotation_data$chr, annotation_data$pos) != paste(ROH_data_pos$chr, ROH_data_pos$pos))
any(paste(annotation_data$chr, annotation_data$pos) != paste(snp_data_pos$chr, snp_data_pos$pos))

# Save all data
save(snp_data, ROH_data, annotation_data, group, file="all_SNP_10KbROH_annotation_DogProjectClareMay2019.Rdata")

