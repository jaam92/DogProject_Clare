
#generate the concatenated and sorted vcfs
#for i in {TenEW,TenIR,TenLB,TenTM,TenBC}
#do 
#/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/vcftools_perl/src/perl/vcf-concat /u/home/c/cmarsden/project-kirk-bigdata/FINAL_VCFS_May2018_INCrels/final_filtered_vcf_justSNPs_6files/"$i"/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_"$i"_jointcalled_chr*.vcf.gz | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/vcftools_perl/src/perl/vcf-sort | bgzip -c > /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ComputeRelatedness/allChroms_"$i"_snpsOnly_sorted.vcf.gz

#for i in {All15AW,All15PG}
#do
#/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/vcftools_perl/src/perl/vcf-concat /u/home/c/cmarsden/project-kirk-bigdata/FINAL_VCFS_May2018_INCrels/final_filtered_vcf_justSNPs_6files/"$i"/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing3_5_v4_mergeGaTKfiltered_varnonvar_"$i"_jointcalled_chr*.vcf.gz | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/vcftools_perl/src/perl/vcf-sort | bgzip -c > /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ComputeRelatedness/allChroms_"$i"_snpsOnly_sorted.vcf.gz


#calculate relatedness using king statistic
#. /u/local/Modules/default/init/modules.sh
#module load vcftools
#for i in {TenEW,TenIR,TenLB,TenTM,TenBC,All15AW,All15PG}
#do
#vcftools --gzvcf allChroms_"$i"_snpsOnly_sorted.vcf.gz  --relatedness2 --out "$i"


#done

#Find First Degree relatives
#remove relatedness with self and drop multiple occurances of same record
#awk '$7 > 0.1 && $1 != $2 {print $0}'  *.relatedness2 | awk '!seen[$1]++' > FirstDegreeRels.txt
