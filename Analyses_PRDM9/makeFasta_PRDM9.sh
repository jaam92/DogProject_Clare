#Load modules
. /u/local/Modules/default/init/modules.sh
module load samtools
module load bcftools

for i in {EW7,AW15,TM3} 
do 

#make fasta
#samtools faidx /u/home/c/cmarsden/project-klohmueldata/clares_data/canids/reference/canFam3/canFam3.fa chr5:63559618-63626468 | bcftools consensus -s "$i" ../master_vcfClare/masterFile_allCanidsN75_ANmt135_Chr5.vcf.gz -o "$i"_PRDM9.fa

#Concatenate
#sed -e 's/chr5:63559618-63626468/'$i'/g' "$i"_PRDM9.fa >> TM3_EW7_AW15_PRDM9.fa

done
