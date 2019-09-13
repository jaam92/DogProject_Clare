library(data.table)
library(dplyr)
library(ggplot2)
library(parallel)
library(stringr)

####IF YOUR NEED TO GENERATE THE INPUT FILE UNCOMMENT
#source(file = "generateRdataFile4Permuation_1Mb_GERP.R")

#Load data
load(file="~/Documents/DogProject_Clare/LocalRscripts/ROHPermutation/all_SNP_1MbROH_annotation_DogProjectClareMay2019.Rdata")
#load(file="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ROHPermutation/all_SNP_1MbROH_annotation_DogProjectClareMay2019.Rdata")
Indivs = read.delim("IndivID.txt")
indNames = Indivs$ID
group = Indivs$Population

# function for calculating the table given the matrices:
calc_contingency_table = function(snpmat, rohmat, annotation, removeFixed=TRUE) {
  
  if (removeFixed) {
    fixedID = apply(snpmat, 1, function(x) all(x==0 | x==-1, na.rm=T) | all(x==2 | x==-1, na.rm=T))    # remove SNPs that are fixed...
    snpmat = snpmat[!fixedID, ] 
    rohmat = rohmat[!fixedID, ]
    annotation = annotation[!fixedID]
  }
  
  rowIDsyn = annotation == "NEU"
  rowIDnonsyn = annotation == "DEL"
  
  temp = (snpmat[rowIDsyn, ])[rohmat[rowIDsyn, ] == 1]
  
  countDerived_SynROH_alleleCopies = sum((temp==1)*1 + (temp==2)*2, na.rm=TRUE)
  countDerived_SynROH_variants = sum((temp==1)*1 + (temp==2)*1, na.rm=TRUE)
  countDerived_SynROH_homozygotes = sum((temp==1)*0 + (temp==2)*1, na.rm=TRUE)
  
  temp = (snpmat[rowIDnonsyn, ])[rohmat[rowIDnonsyn, ] == 1]
  countDerived_NonsynROH_alleleCopies = sum((temp==1)*1 + (temp==2)*2, na.rm=TRUE)
  countDerived_NonsynROH_variants = sum((temp==1)*1 + (temp==2)*1, na.rm=TRUE)
  countDerived_NonsynROH_homozygotes = sum((temp==1)*0 + (temp==2)*1, na.rm=TRUE)
  
  temp = (snpmat[rowIDsyn, ])[rohmat[rowIDsyn, ] == 0]
  countDerived_SynNONROH_alleleCopies = sum((temp==1)*1 + (temp==2)*2, na.rm=TRUE)
  countDerived_SynNONROH_variants = sum((temp==1)*1 + (temp==2)*1, na.rm=TRUE)
  countDerived_SynNONROH_homozygotes = sum((temp==1)*0 + (temp==2)*1, na.rm=TRUE)
  
  temp = (snpmat[rowIDnonsyn, ])[rohmat[rowIDnonsyn, ] == 0]
  countDerived_NonsynNONROH_alleleCopies = sum((temp==1)*1 + (temp==2)*2, na.rm=TRUE)
  countDerived_NonsynNONROH_variants = sum((temp==1)*1 + (temp==2)*1, na.rm=TRUE)
  countDerived_NonsynNONROH_homozygotes = sum((temp==1)*0 + (temp==2)*1, na.rm=TRUE)
  
  table_alleleCopies = matrix(c(countDerived_SynROH_alleleCopies, countDerived_NonsynROH_alleleCopies, countDerived_SynNONROH_alleleCopies, countDerived_NonsynNONROH_alleleCopies), nrow=2, dimnames = list(c("NEU", "DEL"), c("ROH", "nonROH")))
  table_variants = matrix(c(countDerived_SynROH_variants, countDerived_NonsynROH_variants, countDerived_SynNONROH_variants, countDerived_NonsynNONROH_variants), nrow=2, dimnames = list(c("NEU", "DEL"), c("ROH", "nonROH")))
  table_homozygotes = matrix(c(countDerived_SynROH_homozygotes, countDerived_NonsynROH_homozygotes, countDerived_SynNONROH_homozygotes, countDerived_NonsynNONROH_homozygotes), nrow=2, dimnames = list(c("NEU", "DEL"), c("ROH", "nonROH")))
  
  return(list(alleleCopies = table_alleleCopies, variants = table_variants, homozygotes = table_homozygotes))
}

#Calculate empirical table from data
empirical_table = calc_contingency_table(snp_data, ROH_data, annotation_data$annotation, removeFixed = FALSE)

#Compute/look at ratios:
empirical_table$variants[2,1]/empirical_table$variants[1,1]/(empirical_table$variants[2,2]/empirical_table$variants[1,2])

empirical_table$alleleCopies[2,1]/empirical_table$alleleCopies[1,1]/(empirical_table$alleleCopies[2,2]/empirical_table$alleleCopies[1,2])

empirical_table$homozygotes[2,1]/empirical_table$homozygotes[1,1]/(empirical_table$homozygotes[2,2]/empirical_table$homozygotes[1,2])

#Generate tables for each group
empirical_table_group = list()
ratios_group = list()

for (i in  unique(group)){
    empirical_table_group[[i]] = calc_contingency_table(snp_data[,group == i], ROH_data[,group == i], annotation_data$annotation, removeFixed = TRUE)
  
  ratios_group[[i]] = list(variants = empirical_table_group[[i]]$variants[2,1]/empirical_table_group[[i]]$variants[1,1]/(empirical_table_group[[i]]$variants[2,2]/empirical_table_group[[i]]$variants[1,2]),
                             alleleCopies = empirical_table_group[[i]]$alleleCopies[2,1]/empirical_table_group[[i]]$alleleCopies[1,1]/(empirical_table_group[[i]]$alleleCopies[2,2]/empirical_table_group[[i]]$alleleCopies[1,2]),
                             homozygotes = empirical_table_group[[i]]$homozygotes[2,1]/empirical_table_group[[i]]$homozygotes[1,1]/(empirical_table_group[[i]]$homozygotes[2,2]/empirical_table_group[[i]]$homozygotes[1,2]))
}

# Permutations of data (permutation of the SNPs of the ROH matix)
perm_LD = function(matrix) {
  numRow = nrow(matrix)
  shift = sample(1:numRow, size = 1)
  matrixOut = matrix(nrow = numRow, ncol = ncol(matrix))
  matrixOut[1:(numRow-shift+1), ] = matrix[shift:numRow, ]
  matrixOut[(numRow-shift+2):numRow, ] = matrix[1:(shift-1), ]
  return(matrixOut)
}


if(FALSE) {
  jpeg("permutation_1Mb.jpg", width=1200, height=1500/4*5)
  
  par(mfrow=c(5,1))
  image(x=1:length(ROH_data_pos$pos), y=1:ncol(ROH_data), z=as.matrix(ROH_data), col=c("grey", "white", "orange"), xlab="Position", ylab="", yaxt='n')
  axis(2, at = 1:ncol(ROH_data), labels=indNames, las = 1)
  
  image(x=1:length(ROH_data_pos$pos), y=1:ncol(ROH_data), z=matrix(sample(ROH_data, replace = FALSE), ncol=ncol(ROH_data), nrow=nrow(ROH_data)), col=c("grey", "white", "orange"), xlab="Position", ylab="", yaxt='n')
  axis(2, at = 1:ncol(ROH_data), labels=indNames, las = 1)
  
  image(x=1:length(ROH_data_pos$pos), y=1:ncol(ROH_data), z=as.matrix(ROH_data[sample(nrow(ROH_data), replace = FALSE),]), col=c("grey", "white", "orange"), xlab="Position", ylab="", yaxt='n')
  axis(2, at = 1:ncol(ROH_data), labels=indNames, las = 1)
  
  image(x=1:length(ROH_data_pos$pos), y=1:ncol(ROH_data), z=as.matrix(perm_LD(ROH_data)), col=c("grey", "white", "orange"), xlab="Position", ylab="", yaxt='n')
  axis(2, at = 1:ncol(ROH_data), labels=indNames, las = 1)
  
  snp_data_ann = snp_data
  snp_data_ann[annotation_data$annotation=="DEL",] = snp_data_ann[annotation_data$annotation=="DEL",] + 4
  
  image(x=1:length(snp_data_pos$pos), y=1:ncol(snp_data), z=as.matrix(snp_data_ann), col=c("white", "white", "lightblue", "blue", "white", "white", "coral", "red"), xlab="Position", ylab="", yaxt='n')
  axis(2, at = 1:ncol(ROH_data), labels=indNames, las = 1)
  
  dev.off()
}

#Permutation test full table across countries (permutating columns):
#library(parallel)
pval_group = list()
X_squared_group = list()

nrep = 10000
for (i in unique(group)){
  print(i)
    ID = group == i
  empirical_table = empirical_table_group[[i]]
  permutation_tables_Columns = mclapply(1:nrep, function(x) calc_contingency_table(snp_data[,ID], (ROH_data[,ID])[sample(nrow(ROH_data[,ID]), replace = FALSE),], annotation_data$annotation, removeFixed = TRUE), mc.cores = 4)
  X_squared_nulldistr_Columns = data.frame(t(sapply(permutation_tables_Columns, function(x) sapply(x, function(y) chisq.test(y)$statistic))))
  
  X_squared_group[[i]] = list(variants = X_squared_nulldistr_Columns$variants.X.squared,
                              alleleCopies = X_squared_nulldistr_Columns$alleleCopies.X.squared,
                              homozygotes = X_squared_nulldistr_Columns$homozygotes.X.squared)
  
  pval_group[[i]] = list(variants = sum(X_squared_nulldistr_Columns$variants.X.squared >= chisq.test(empirical_table$variants)$statistic)/nrep,
                           alleleCopies = sum(X_squared_nulldistr_Columns$alleleCopies.X.squared >= chisq.test(empirical_table$alleleCopies)$statistic)/nrep,
                           homozygotes = sum(X_squared_nulldistr_Columns$homozygotes.X.squared >= chisq.test(empirical_table$homozygotes)$statistic)/nrep)
}


save(file="permutation_results_allCanids_1MbROH_GERP.Rdata", pval_group, X_squared_group, empirical_table_group, ratios_group)



#Plots

load(file="permutation_results_allCanids_1MbROH_GERP.Rdata")

ratios_group_df = rbind(data.frame(type = "variants", ratio = sapply(ratios_group, function(x) x$variants)), 
                          data.frame(type = "alleleCopies", ratio = sapply(ratios_group, function(x) x$alleleCopies)), 
                          data.frame(type = "homozygotes", ratio = sapply(ratios_group, function(x) x$homozygotes)))

ratios_group_df$group = rownames(ratios_group_df)[1:7]

pval_group_df = rbind(data.frame(type = "variants", pvalue = sapply(pval_group, function(x) x$variants)), 
                        data.frame(type = "alleleCopies", pvalue = sapply(pval_group, function(x) x$alleleCopies)), 
                        data.frame(type = "homozygotes", pvalue = sapply(pval_group, function(x) x$homozygotes)))

pval_group_df$group = rownames(pval_group_df)[1:7]

ratios_group_df$pvalue = pval_group_df$pvalue

#ggplot(data = ratios_group_df) + geom_point(aes(group, ratio, col=type, cex=pvalue<0.05)) + geom_line(aes(group, ratio)) + theme_bw()

write.csv(ratios_group_df, row.names = FALSE, file = "ratios_pvalues_pergroup_1MbROH_GERP.csv")
