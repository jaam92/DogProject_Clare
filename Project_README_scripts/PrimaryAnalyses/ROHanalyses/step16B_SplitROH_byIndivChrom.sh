

#Split by chrom and indiv
awk '{print $2"\t"$3 >  "FinalROH_"$11"_"$1"_1MbRanges.txt"}' ../FinalROHQCd_min1Mb_allIndivs_Dec2018_DogProjClare.LROH
