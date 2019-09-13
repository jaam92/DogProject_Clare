
#Load modules
. /u/local/Modules/default/init/modules.sh
module load plink/1.90b3.45
#module load python/3.7.0

#vcf to plink file
#plink --dog --vcf ../master_vcfClare/masterFile_allChroms.vcf --make-bed --out masterFile_allChroms

#generate rsIDs because this is the identifier for LD pruning and our vcf does not have rsIDs
#python3 genRSid.py > random_rsIDs.txt

#rename the snps
#paste masterFile_allChroms.bim random_rsIDs.txt | awk '{print $1"\t"$7"\t"$3"\t"$4"\t"$5"\t"$6}' > masterFile_allChroms_v2.bim
#rm masterFile_allChroms.bim
#mv masterFile_allChroms_v2.bim masterFile_allChroms.bim 

#ld prune the plink file
#Pruning for ADMIXTURE
#plink --dog --bfile masterFile_allChroms --indep-pairwise 50 10 0.1

#Pruning for smartPCA
#plink --dog --bfile masterFile_allChroms --maf 0.05 --indep-pairwise 50 5 0.8

#make pruned plink file
#plink --dog --bfile masterFile_allChroms --extract plink.prune.in --make-bed --out prunedData_masterFile_allChroms_ADMIXTURE

plink --dog --bfile masterFile_allChroms --extract plink.prune.in --make-bed --out prunedData_masterFile_allChroms_smartPCA

