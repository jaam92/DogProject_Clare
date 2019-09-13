


#Generate the ROH file but separate by species and/or breed
. /u/local/Modules/default/init/modules.sh
module load vcftools

for f in {1..38}
do

vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allCanidsN75_ANmt135_Chr"$f".vcf.gz --LROH --chr chr"$f" --out /u/flashscratch/j/jmooney3/vcfToolsROH/jackieVersion_Chr"$f"

done

