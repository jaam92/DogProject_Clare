#Load modules
. /u/local/Modules/default/init/modules.sh
module load samtools
module load bcftools

for i in {EW7,AW15,TM3} 
do 

#generate the fasta file
#samtools faidx /u/home/c/cmarsden/project-klohmueldata/clares_data/canids/reference/canFam3/canFam3.fa chr27:38468698-38472776 | bcftools consensus -s "$i" ../master_vcfClare/masterFile_allCanidsN75_ANmt135_Chr27.vcf.gz -o "$i"_GAPDH.fa


#Concatenate
#sed -e 's/chr27:38468698-38472776/'$i'/g' "$i"_GAPDH.fa >> TM3_EW7_AW15_GAPDH.fa

done
