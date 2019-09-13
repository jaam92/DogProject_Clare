
##qsub -cwd -V -N tpedformat -l highp,time=02:30:00,h_data=4G -t 1-38:1 -M eplau -m bea vcf_to_tped.sh

##Deleted tped files after plink binaries were made

#set chroms as num of tasks
CHROM=${SGE_TASK_ID}

#load vcftools
. /u/local/Modules/default/init/modules.sh
module load vcftools

#Specify directories
out_dir='/u/flashscratch/j/jmooney3/DogsROH/tped_Merged'
in_dir='/u/flashscratch/j/jmooney3/DogsROH/masterVCFs'

#Specify files
merged_plinkout=${out_dir}/'allSpecies_MasterFile_varnonvar_jointcalled_biallelic_chr'${CHROM}''
merged_vcfin=${in_dir}/'allSpecies_MasterFile_varnonvar_jointcalled_biallelic_chr'${CHROM}'.vcf.gz'

#Convert vcf to tped
vcftools --gzvcf ${merged_vcfin} --plink-tped --out ${merged_plinkout}

sleep 200
