##qsub -cwd -V -N genMasterFile -l highp,time=18:00:00,h_data=12G -t 1-38:1 -M eplau -m a createAllSitesVCFs_N6.sh

in_dir='/u/home/c/cmarsden/project-kirk-bigdata/FINAL_VCFS_May2018_INCrels/final_master_allsites_6files_vcf'
out_dir='/u/flashscratch/j/jmooney3/AllSitesFiles/SubsetN6'
indivFile_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles'
gatk_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/GATK/gatk_35'
ref_genome='/u/home/c/cmarsden/project-klohmueldata/clares_data/canids/reference/canFam3/canFam3.fa'
MEM=4G

#set chroms as num of tasks
CHROM=${SGE_TASK_ID}

#infile
vcfin=${in_dir}/'MasterEVERYTHING_7gps_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2Ten_3All15_5_v4_mergeGaTKfiltered_varnonvar_TenBCLBTMIREWAll15PGAW_jointcalled_chr'${CHROM}'.vcf.gz'

#outfiles
EW_PostFilter=${out_dir}/'EW_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
IR_PostFilter=${out_dir}/'IR_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
AW_PostFilter=${out_dir}/'AW_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
TM_PostFilter=${out_dir}/'TM_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
LB_PostFilter=${out_dir}/'LB_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
PG_PostFilter=${out_dir}/'PG_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'
BC_PostFilter=${out_dir}/'BC_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz'

#Individuals
keep_EW=${indivFile_dir}/'EW_indivList_downsample_N6.txt'
keep_IR=${indivFile_dir}/'IR_indivList_downsample_N6.txt'
keep_AW=${indivFile_dir}/'AW_indivList_downsample_N6.txt'
keep_TM=${indivFile_dir}/'TM_indivList_downsample_N6.txt'
keep_LB=${indivFile_dir}/'LB_indivList_downsample_N6.txt'
keep_PG=${indivFile_dir}/'PG_indivList_downsample_N6.txt'
keep_BC=${indivFile_dir}/'BC_indivList_downsample_N6.txt'

#submit script
. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77

########Subset Indivs and subset to sites with 100% coverage########

java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
	-T SelectVariants \
	--variant ${vcfin} \
	--sample_file ${keep_EW} \
	-select "AN > 11" \
	-o ${EW_PostFilter}

java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
        -T SelectVariants \
        --variant ${vcfin} \
        --sample_file ${keep_IR} \
        -select "AN > 11" \
        -o ${IR_PostFilter}

java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
        -T SelectVariants \
        --variant ${vcfin} \
        --sample_file ${keep_AW} \
        -select "AN > 11" \
        -o ${AW_PostFilter}

java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
        -T SelectVariants \
        --variant ${vcfin} \
        --sample_file ${keep_TM} \
        -select "AN > 11" \
        -o ${TM_PostFilter}

java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
        -T SelectVariants \
        --variant ${vcfin} \
        --sample_file ${keep_LB} \
        -select "AN > 11" \
        -o ${LB_PostFilter}

java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
        -T SelectVariants \
        --variant ${vcfin} \
        --sample_file ${keep_PG} \
        -select "AN > 11" \
        -o ${PG_PostFilter}

java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
        -T SelectVariants \
        --variant ${vcfin} \
        --sample_file ${keep_BC} \
        -select "AN > 11" \
        -o ${BC_PostFilter}
