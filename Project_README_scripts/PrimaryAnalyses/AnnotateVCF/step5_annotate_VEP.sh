##qsub -cwd -V -N VEP -l highp,time=05:00:00,h_data=1G -t 1-38:1 -M eplau -m a annotate_VEP.sh

#set chroms as num of tasks
CHROM=${SGE_TASK_ID}
#CHROM=38

##Directories
VEP_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/ensembl-vep-release-94'
in_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare'
out_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/VEP'

#Input and output
vcf_in=${in_dir}/'masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.vcf.gz'
vcf_out=${out_dir}/'annotatedVEPandSIFT_masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.txt'
stats_fname=${out_dir}/'statFile_masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.html'

#submit script
. /u/local/Modules/default/init/modules.sh
module load perl

perl ${VEP_dir}/vep --tab -database -i ${vcf_in} -o ${vcf_out} --stats_file ${stats_fname} --sift s --species canis_familiaris --allow_non_variant --symbol --canonical --force_overwrite --fields "Location,Allele,Consequence,IMPACT,BIOTYPE,CANONICAL,SIFT,SYMBOL,Gene"

