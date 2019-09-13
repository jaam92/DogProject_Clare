#Load modules
. /u/local/Modules/default/init/modules.sh
module load samtools
module load bcftools

for i in {EW7,AW15,TM3,BC6} 
do
#generate fasta from vcf for single individual
#samtools faidx /u/home/c/cmarsden/project-klohmueldata/clares_data/canids/reference/canFam3/canFam3.fa chr6:37409924-37538007 | bcftools consensus -s "$i" ../master_vcfClare/masterFile_allCanidsN75_ANmt135_Chr6.vcf.gz -o "$i"_CREBBP.fa

#Concatenate
#sed -e 's/chr6:37409924-37538007/'$i'/g' "$i"_CREBBP.fa >> TM3_EW7_AW15_BC6_CREBBP.fa 

done
