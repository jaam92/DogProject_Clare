#Takes about 3 minutes to run


for f in {1..38} 
do 

#convert tped to binaries
plink --tfile tped_Merged/allSpecies_MasterFile_varnonvar_jointcalled_biallelic_chr"$f" --make-bed --dog --out mergedPLINK/allSpecies_MasterFile_varnonvar_jointcalled_biallelic_chr"$f"

done
