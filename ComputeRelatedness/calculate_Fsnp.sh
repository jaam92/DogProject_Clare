


for i in {TenEW,TenIR,TenLB,TenTM,TenBC,All15AW,All15PG} 
do 

#Calculate Fsnp using common sites (MAF >10%)
#. /u/local/Modules/default/init/modules.sh
#module load vcftools
#vcftools --gzvcf allChroms_"$i"_snpsOnly_sorted.vcf.gz --maf 0.1 --het --out Fsnp_allChroms_"$i"_snpsOnly_sorted

done

#Concatenate and remove repeats of the header line
#cat *_sorted.het | awk '!seen[$1] ++' > allPops_Fsnp_n81.het
