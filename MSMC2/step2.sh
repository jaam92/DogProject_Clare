##qsub -cwd -V -N createMask -l time=03:00:00,h_data=8G -t 1-38:1 -M eplau -m ea step2.sh
in_dir='/u/flashscratch/j/jmooney3/MSMC2Inputs'
out_dir='/u/flashscratch/j/jmooney3/MSMC2Inputs'
#out_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/MSMC2'
software_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msmc2-master/msmc-tools-master'

#set chroms as num of tasks
CHROM=${SGE_TASK_ID}

#infile
EW2_FIN=${in_dir}/'EW2_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
EW3_FIN=${in_dir}/'EW3_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
EW9_FIN=${in_dir}/'EW9_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
EW10_FIN=${in_dir}/'EW10_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
AW11_FIN=${in_dir}/'AW11_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
AW20_FIN=${in_dir}/'AW20_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
IR5_FIN=${in_dir}/'IR5_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
IR10_FIN=${in_dir}/'IR10_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'

#mask outfiles 
EW2_mask=${out_dir}/'EW2_chr'${CHROM}'Mask.bed.gz'
EW3_mask=${out_dir}/'EW3_chr'${CHROM}'Mask.bed.gz'
EW9_mask=${out_dir}/'EW9_chr'${CHROM}'Mask.bed.gz'
EW10_mask=${out_dir}/'EW10_chr'${CHROM}'Mask.bed.gz'
AW11_mask=${out_dir}/'AW11_chr'${CHROM}'Mask.bed.gz'
AW20_mask=${out_dir}/'AW20_chr'${CHROM}'Mask.bed.gz'
IR5_mask=${out_dir}/'IR5_chr'${CHROM}'Mask.bed.gz'
IR10_mask=${out_dir}/'IR10_chr'${CHROM}'Mask.bed.gz'

#outfiles
EW2_FOUT=${out_dir}/'EW2_chr'${CHROM}'_output.vcf'
EW3_FOUT=${out_dir}/'EW3_chr'${CHROM}'_output.vcf'
EW9_FOUT=${out_dir}/'EW9_chr'${CHROM}'_output.vcf'
EW10_FOUT=${out_dir}/'EW10_chr'${CHROM}'_output.vcf'
AW11_FOUT=${out_dir}/'AW11_chr'${CHROM}'_output.vcf'
AW20_FOUT=${out_dir}/'AW20_chr'${CHROM}'_output.vcf'
IR5_FOUT=${out_dir}/'IR5_chr'${CHROM}'_output.vcf'
IR10_FOUT=${out_dir}/'IR10_chr'${CHROM}'_output.vcf'

#outfiles zipped
EW2_FOUT_ZIPPED=${out_dir}/'EW2_chr'${CHROM}'_output.vcf.gz'
EW3_FOUT_ZIPPED=${out_dir}/'EW3_chr'${CHROM}'_output.vcf.gz'
EW9_FOUT_ZIPPED=${out_dir}/'EW9_chr'${CHROM}'_output.vcf.gz'
EW10_FOUT_ZIPPED=${out_dir}/'EW10_chr'${CHROM}'_output.vcf.gz'
AW11_FOUT_ZIPPED=${out_dir}/'AW11_chr'${CHROM}'_output.vcf.gz'
AW20_FOUT_ZIPPED=${out_dir}/'AW20_chr'${CHROM}'_output.vcf.gz'
IR5_FOUT_ZIPPED=${out_dir}/'IR5_chr'${CHROM}'_output.vcf.gz'
IR10_FOUT_ZIPPED=${out_dir}/'IR10_chr'${CHROM}'_output.vcf.gz'

#Load python
. /u/local/Modules/default/init/modules.sh
module load python/2.7

#Create mask and new vcf
zcat ${EW2_FIN} | python ${software_dir}/vcfAllSiteParser.py chr${CHROM} ${EW2_mask} > ${EW2_FOUT}
zcat ${EW3_FIN} | python ${software_dir}/vcfAllSiteParser.py chr${CHROM} ${EW3_mask} > ${EW3_FOUT}
zcat ${EW9_FIN} | python ${software_dir}/vcfAllSiteParser.py chr${CHROM} ${EW9_mask} > ${EW9_FOUT}
zcat ${EW10_FIN} | python ${software_dir}/vcfAllSiteParser.py chr${CHROM} ${EW10_mask} > ${EW10_FOUT}
zcat ${AW11_FIN} | python ${software_dir}/vcfAllSiteParser.py chr${CHROM} ${AW11_mask} > ${AW11_FOUT}
zcat ${AW20_FIN} | python ${software_dir}/vcfAllSiteParser.py chr${CHROM} ${AW20_mask} > ${AW20_FOUT}
zcat ${IR5_FIN} | python ${software_dir}/vcfAllSiteParser.py chr${CHROM} ${IR5_mask} > ${IR5_FOUT}
zcat ${IR10_FIN} | python ${software_dir}/vcfAllSiteParser.py chr${CHROM} ${IR10_mask} > ${IR10_FOUT}

#bgzip and tabix index
bgzip -c ${EW2_FOUT} > ${EW2_FOUT_ZIPPED}
bgzip -c ${EW3_FOUT} > ${EW3_FOUT_ZIPPED}
bgzip -c ${EW9_FOUT} > ${EW9_FOUT_ZIPPED}
bgzip -c ${EW10_FOUT} > ${EW10_FOUT_ZIPPED}
bgzip -c ${AW11_FOUT} > ${AW11_FOUT_ZIPPED}
bgzip -c ${AW20_FOUT} > ${AW20_FOUT_ZIPPED}
bgzip -c ${IR5_FOUT} > ${IR5_FOUT_ZIPPED}
bgzip -c ${IR10_FOUT} > ${IR10_FOUT_ZIPPED}

tabix -p vcf ${EW2_FOUT_ZIPPED}
tabix -p vcf ${EW3_FOUT_ZIPPED}
tabix -p vcf ${EW9_FOUT_ZIPPED}
tabix -p vcf ${EW10_FOUT_ZIPPED}
tabix -p vcf ${AW11_FOUT_ZIPPED}
tabix -p vcf ${AW20_FOUT_ZIPPED}
tabix -p vcf ${IR5_FOUT_ZIPPED}
tabix -p vcf ${IR10_FOUT_ZIPPED}

