##This Rscript will call Nandita's H12 script and run it on the file set

#load packages
library(reticulate)
library(dplyr)
library(readr)
library(stringr)
library(magrittr)
library(glue)

setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/ComputeG12/RunAnalyses")

#Step 1 is to generate the haplotype files
#Use default from nandita
fnames = dir(path = "/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/ComputeG12", pattern = "\\inputG12haps.txt$")

for (i in seq_along(fnames)){
  inputFile = fnames[i]
  population = gsub('_n.*','',inputFile)
  numSamps = parse_number(str_extract(inputFile, "_n[^_]*"))
  chromo = parse_number(str_extract(inputFile, "Chr[^_]*"))
  
  pythonInputs = glue('
                      /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/ComputeG12/{inputFile} {numSamps} -o output_{inputFile} -w 200 -j 50 -d 0
                      ')
  
  pathScript = glue('python /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/SelectionHapStats-master/scripts/H12_H2H1.py')
  
  command1 = glue('{pathScript} {pythonInputs}')
  
  system(command1)
  
  message1 = sprintf("running input file %s", inputFile)
  print(message1)
  
}


#Step 2 is to find the selection peaks in the haplotype files
#Use the default from nandita
hapFiles = dir(path = "/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/ComputeG12/RunAnalyses/", pattern = "^output\\w*\\.txt$")

for (j in seq_along(hapFiles)){
  inputFile2 = hapFiles[j]
  
  pythonInputs2 = glue('
                      /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/ComputeG12/RunAnalyses/{inputFile2} -o Peaks_{inputFile2} -t 0.02 
                      ')
  
  pathScript2 = glue('python /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/SelectionHapStats-master/scripts/H12peakFinder.py')
  
  command2 = glue('{pathScript2} {pythonInputs2}')
  
  system(command2)
  
  message2 = sprintf("running input file %s", inputFile2)
  print(message2)
}

