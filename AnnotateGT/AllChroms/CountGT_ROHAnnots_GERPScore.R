####This script will take in annotation data across all chromosome for each individual and count up the types of annotations (with three different counting methods), how many annotations fall within or outside of ROH of a given length, compute OR of NS to SYN variants within vs outside of ROH for each individual and across populations

####Input file: annotation file for each individual with data from all chroms
####Output file: File with above information for all individuals and a file with OR for different counting methods split by population 

#Load Libraries
library(dplyr)
library(data.table)
library(mgsub)
library(reshape2)
library(tidyr)

#Get input files
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT/AllChroms")
filenames = list.files(pattern = glob2rx("annotatedGTwithVEP_*_allChroms.txt"))

#Empty data frame to fill with data
allIndivDF = data.frame()

#Loop thru all individuals files
for (i in 1:length(filenames)){
  
  #Read files in 
  annots = read.delim(file = filenames[i])
  
  #Grab individual ID
  iid = mgsub_dict(filenames[i], conversions = list("annotatedGTwithVEP_" = "", "_allChroms.txt"="")) 
  
####Grab counts ROH of at least 10Kb
  CountDF_withROH = annots %>%
    filter(GT != "AncHom" & GT != "Missing") %>% #Remove AncHom and Missing
    group_by(GT,ANNOT,withinROH_10Kb) %>%
    count() %>% #Count Variants (DERHOM =1 , HET = 1)
    mutate(CountAlleles = ifelse(GT == "DerHom", as.numeric(n)*2, as.numeric(n)),#Count Alleles (DERHOM =2 , HET = 1)
           CountDerHom = ifelse(GT == "DerHom", as.numeric(n), as.numeric(0))) %>% #Count DerHOM (DERHOM =1 )
    group_by(ANNOT,withinROH_10Kb) %>%
    summarise_at(.vars= vars(n, CountAlleles, CountDerHom), .funs = sum) %>%
    plyr::rename(c("n"="CountVariants"))

####Grab counts ROH of at least 1Mb
  CountDF_withROH_1Mb = annots %>%
    filter(GT != "AncHom" & GT != "Missing") %>% 
    group_by(GT,ANNOT,withinROH_1Mb) %>%
    count() %>% 
    mutate(CountAlleles = ifelse(GT == "DerHom", as.numeric(n)*2, as.numeric(n)),
           CountDerHom = ifelse(GT == "DerHom", as.numeric(n), as.numeric(0))) %>%
    group_by(ANNOT,withinROH_1Mb) %>%
    summarise_at(.vars= vars(n, CountAlleles, CountDerHom), .funs = sum) %>%
    plyr::rename(c("n"="CountVariants"))

####Grab counts within and outside ROH of at least 10Kb
  ROHStats_10Kb = melt(CountDF_withROH, id.vars=c("withinROH_10Kb", "ANNOT")) %>% 
    arrange(ANNOT) %>%
    mutate(withinROH_10Kb = ifelse(withinROH_10Kb == 1, "ROH", "nonROH"),
           NewColName = paste(ANNOT,"_",variable,"_",withinROH_10Kb,sep = "")) %>%
    select(NewColName,value) %>%
    spread(NewColName,value)
####Grab counts within and outside ROH of at least 1Mb  
  ROHStats_1Mb = melt(CountDF_withROH_1Mb, id.vars=c("withinROH_1Mb", "ANNOT")) %>% 
    arrange(ANNOT) %>%
    mutate(withinROH_1Mb = ifelse(withinROH_1Mb == 1, "ROH1Mb", "nonROH1Mb"),
           NewColName = paste(ANNOT,"_",variable,"_",withinROH_1Mb,sep = "")) %>%
    select(NewColName,value) %>%
    spread(NewColName,value)

###Counts of just annotations
  CountDF_withJustAnnots = CountDF_withROH %>% 
    group_by(ANNOT) %>% 
    summarise_at(.vars= vars(CountVariants, CountAlleles, CountDerHom), .funs = sum) 
  
####Counts of variants by Impact
  ImpactCounts = annots %>%
    group_by(IMPACT) %>%
    count() %>%
    spread(IMPACT,n)
  
####Count up the total number of annotated variants and how they are distributed
  Calls = annots %>%
    group_by(GT) %>%
    count() %>%
    spread(GT,n) %>%
    mutate(LineCount = dim(annots)[1])
  
####Counts GERP Scores
  GERPScore = annots %>%
    filter(GERPScore >= 4 | GERPScore <= 2, GT != "AncHom", GT != "Missing") %>%
    filter(GERPScore <= 2 | GERPScore >= 4 & ANNOT != "SY") %>% #Remove SY variation GERP says is deleterious
    mutate(GERPScore = ifelse(GERPScore >= 4, "PutDel", "PutNeu")) %>%
    group_by(GT,GERPScore) %>%
    count() %>%
    mutate(CountAlleles = ifelse(GT == "DerHom", as.numeric(n)*2, as.numeric(n)),
           CountDerHom = ifelse(GT == "DerHom", as.numeric(n), as.numeric(0))) %>%
    group_by(GERPScore) %>%
    summarise_at(.vars= vars(n, CountAlleles, CountDerHom), .funs = sum) %>%
    plyr::rename(c("n"="CountVariants")) %>%
    melt(id.vars=c("GERPScore")) %>% #convert table from wide to long
    arrange(GERPScore) %>% #order by annotation
    mutate(NewColName = paste(GERPScore,"_",variable, sep = "")) %>%
    select(NewColName,value) %>% #only keep new col and values
    spread(NewColName,value) 

####Put all the counts together and add odds-ratios within vs outside ROH  
  indivLongDF = melt(CountDF_withJustAnnots, id.vars=c("ANNOT")) %>% #convert table from wide to long
    arrange(ANNOT) %>% #order by annotation
    mutate(NewColName = paste(ANNOT,"_",variable, sep = "")) %>%
    select(NewColName,value) %>% #only keep new col and values
    spread(NewColName,value) %>% #make data one row
    cbind(Calls,ImpactCounts,GERPScore,ROHStats_10Kb,ROHStats_1Mb) %>% #add on total calls and impact counts
    mutate(ID = iid, 
           Population = substr(iid,1,2),
           OR_AlleleCopies = (SY_CountAlleles_nonROH*NS_CountAlleles_ROH)/(SY_CountAlleles_ROH*NS_CountAlleles_nonROH),
           OR_Variants = (SY_CountVariants_nonROH*NS_CountVariants_ROH)/(SY_CountVariants_ROH*NS_CountVariants_nonROH),
           OR_DerHom = (SY_CountDerHom_nonROH*NS_CountDerHom_ROH)/(SY_CountDerHom_ROH*NS_CountDerHom_nonROH),
           OR_AlleleCopies_ROH1Mb = (SY_CountAlleles_nonROH1Mb*NS_CountAlleles_ROH1Mb)/(SY_CountAlleles_ROH1Mb*NS_CountAlleles_nonROH1Mb),
           OR_Variants_ROH1Mb = (SY_CountVariants_nonROH1Mb*NS_CountVariants_ROH1Mb)/(SY_CountVariants_ROH1Mb*NS_CountVariants_nonROH1Mb),
           OR_DerHom_ROH1Mb = (SY_CountDerHom_nonROH1Mb*NS_CountDerHom_ROH1Mb)/(SY_CountDerHom_ROH1Mb*NS_CountDerHom_nonROH1Mb)) %>% #add Individual ID, Population, and OR for counting methods and different types of ROH
    select(ID, Population, everything())#move individual id and population to front
  
  #Dataframe with all individuals (use bind rows so columns w/out values filled as na)
  allIndivDF = bind_rows(allIndivDF, indivLongDF)
  
}

#Write output file
write.table(allIndivDF, file="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT/AllChroms/GTAnnotationCountResults_May2019_DogProjClare.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Calculate OR for entire population
PerPopulationOR = allIndivDF %>%
  group_by(Population) %>%
  summarise(OR_AlleleCopies =(sum(SY_CountAlleles_nonROH)*sum(NS_CountAlleles_ROH))/(sum(SY_CountAlleles_ROH)*sum(NS_CountAlleles_nonROH)),
            OR_Variants = (sum(SY_CountVariants_nonROH)*sum(NS_CountVariants_ROH))/(sum(SY_CountVariants_ROH)*sum(NS_CountVariants_nonROH)),
            OR_DerHom = (sum(SY_CountDerHom_nonROH)*sum(NS_CountDerHom_ROH))/(sum(SY_CountDerHom_ROH)*sum(NS_CountDerHom_nonROH)),
            OR_AlleleCopies_1Mb =(sum(SY_CountAlleles_nonROH1Mb)*sum(NS_CountAlleles_ROH1Mb))/(sum(SY_CountAlleles_ROH1Mb)*sum(NS_CountAlleles_nonROH1Mb)),
            OR_Variants_1Mb = (sum(SY_CountVariants_nonROH1Mb)*sum(NS_CountVariants_ROH1Mb))/(sum(SY_CountVariants_ROH1Mb)*sum(NS_CountVariants_nonROH1Mb)),
            OR_DerHom_1Mb = (sum(SY_CountDerHom_nonROH1Mb)*sum(NS_CountDerHom_ROH1Mb))/(sum(SY_CountDerHom_ROH1Mb)*sum(NS_CountDerHom_nonROH1Mb))) %>%
  mutate(Population = mgsub_dict(Population, conversions =  list("BC" = "Border Collie", "LB" = "Labrador Retriever", "PG" = "Pug", "TM" = "Tibetan Mastiff", "AW" = "Arctic Wolf", "EW" = "Ethiopian Wolf", "IR" = "Isle Royale"))) #Add population name back

#Write output file
write.table(PerPopulationOR, file="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/AnnotateGT/AllChroms/PerPopulationOR_May2019_DogProjClare.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


