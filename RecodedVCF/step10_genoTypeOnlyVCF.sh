source /u/local/Modules/default/init/modules.sh
module load bcftools

for f in {1..38} 
do 

#create a vcf with (CHROM,POS,REF,ALT,GT)
#remove all the extraneous info in the indiv column

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allCanidsN75_ANmt135_Chr"$f".vcf.gz -o /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/GenotypeOnly/reformatted_masterFile_allCanidsN75_ANmt135_Chr"$f".vcf

done 


#remove weird space between the hash and CHROM
#sed -i 's/# /#/g' reformatted_masterFile_allCanidsN75_ANmt135_Chr*.vcf
