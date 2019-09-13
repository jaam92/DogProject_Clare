##qsub -cwd -V -N genMSMCFile -l time=01:00:00,h_data=1G -t 1-38:1 -M eplau -m a step1_createSingleIndivVCF.sh
in_dir='/u/flashscratch/j/jmooney3/MSMC2Inputs'
out_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/MSMC2'
software_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/msmc2-master/msmc-tools-master'

#set chroms as num of tasks
CHROM=${SGE_TASK_ID}

#mask outfiles 
EW2_mask=${in_dir}/'EW2_chr'${CHROM}'Mask.bed.gz'
EW3_mask=${in_dir}/'EW3_chr'${CHROM}'Mask.bed.gz'
EW9_mask=${in_dir}/'EW9_chr'${CHROM}'Mask.bed.gz'
EW10_mask=${in_dir}/'EW10_chr'${CHROM}'Mask.bed.gz'
AW11_mask=${in_dir}/'AW11_chr'${CHROM}'Mask.bed.gz'
AW20_mask=${in_dir}/'AW20_chr'${CHROM}'Mask.bed.gz'
IR5_mask=${in_dir}/'IR5_chr'${CHROM}'Mask.bed.gz'
IR10_mask=${in_dir}/'IR10_chr'${CHROM}'Mask.bed.gz'

#infiles
EW2_FIN=${in_dir}/'EW2_chr'${CHROM}'_output.vcf.gz'
EW3_FIN=${in_dir}/'EW3_chr'${CHROM}'_output.vcf.gz'
EW9_FIN=${in_dir}/'EW9_chr'${CHROM}'_output.vcf.gz'
EW10_FIN=${in_dir}/'EW10_chr'${CHROM}'_output.vcf.gz'
AW11_FIN=${in_dir}/'AW11_chr'${CHROM}'_output.vcf.gz'
AW20_FIN=${in_dir}/'AW20_chr'${CHROM}'_output.vcf.gz'
IR5_FIN=${in_dir}/'IR5_chr'${CHROM}'_output.vcf.gz'
IR10_FIN=${in_dir}/'IR10_chr'${CHROM}'_output.vcf.gz'

#outfiles
EW2_FOUT=${out_dir}/'EW2_chr'${CHROM}'_msmcInput.txt'
EW3_FOUT=${out_dir}/'EW3_chr'${CHROM}'_msmcInput.txt'
EW9_FOUT=${out_dir}/'EW9_chr'${CHROM}'_msmcInput.txt'
EW10_FOUT=${out_dir}/'EW10_chr'${CHROM}'_msmcInput.txt'
AW11_FOUT=${out_dir}/'AW11_chr'${CHROM}'_msmcInput.txt'
AW20_FOUT=${out_dir}/'AW20_chr'${CHROM}'_msmcInput.txt'
IR5_FOUT=${out_dir}/'IR5_chr'${CHROM}'_msmcInput.txt'
IR10_FOUT=${out_dir}/'IR10_chr'${CHROM}'_msmcInput.txt'

#Load python
. /u/local/Modules/default/init/modules.sh
module load python/3.4

#Generate msmc input
python3 ${software_dir}/generate_multihetsep.py --mask ${EW2_mask} ${EW2_FIN} > ${EW2_FOUT}
python3 ${software_dir}/generate_multihetsep.py --mask ${EW3_mask} ${EW3_FIN} > ${EW3_FOUT}
python3 ${software_dir}/generate_multihetsep.py --mask ${EW9_mask} ${EW9_FIN} > ${EW9_FOUT}
python3 ${software_dir}/generate_multihetsep.py --mask ${EW10_mask} ${EW10_FIN} > ${EW10_FOUT}
python3 ${software_dir}/generate_multihetsep.py --mask ${AW11_mask} ${AW11_FIN} > ${AW11_FOUT}
python3 ${software_dir}/generate_multihetsep.py --mask ${AW20_mask} ${AW20_FIN} > ${AW20_FOUT}
python3 ${software_dir}/generate_multihetsep.py --mask ${IR5_mask} ${IR5_FIN} > ${IR5_FOUT}
python3 ${software_dir}/generate_multihetsep.py --mask ${IR10_mask} ${IR10_FIN} > ${IR10_FOUT}


sleep 200
