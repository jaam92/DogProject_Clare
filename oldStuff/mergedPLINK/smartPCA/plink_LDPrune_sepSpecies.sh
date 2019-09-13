#LD prune for each species
#for i in {AW,EW,IR,GS,PG,TM}; do plink --bfile allDogsMerged_smartPCA --keep "$i"_indivs.txt --maf 0.05 --indep-pairwise 50 5 0.5 --dog --out sepBySpecies/"$i"_keepSites_LDPrune;done


#Make new bed file separated by species and LD pruned
#for i in {AW,EW,IR,GS,PG,TM}; do plink --bfile allDogsMerged_smartPCA --keep "$i"_indivs.txt --extract sepBySpecies/"$i"_keepSites_LDPrune.prune.in --make-bed --dog --out sepBySpecies/"$i"_smartPCA_LDpruned;done
