##qsub -cwd -V -N genMasterFile -l highp,time=02:30:00,h_data=10G -t 1-38:1 -M eplau -m ea createMasterVCFs.sh

in_dir='/u/home/c/cmarsden/project-kirk-bigdata/FINAL_VCFS_May2018_INCrels/final_master_allsites_6files_vcf'
out_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare'
gatk_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/GATK/gatk_35'
ref_genome='/u/home/c/cmarsden/project-klohmueldata/clares_data/canids/reference/canFam3/canFam3.fa'
MEM=4G

#set chroms as num of tasks
CHROM=${SGE_TASK_ID}

#infile contains all sites from all breeds that passed filtering
vcfin=${in_dir}/'MasterEVERYTHING_7gps_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2Ten_3All15_5_v4_mergeGaTKfiltered_varnonvar_TenBCLBTMIREWAll15PGAW_jointcalled_chr'${CHROM}'.vcf.gz'

#outfiles are sites that are biallelic in at least one individual
Canids_PostFilter=${out_dir}/'masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.vcf'

#Individuals
keep_Canids='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/Dogs2Keep.txt'

#submit script
. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77

########Subset Indivs and subset to sites with 90% coverage########

java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
	-T SelectVariants \
	--variant ${vcfin} \
	--sample_file ${keep_Canids} \
	-selectType SNP \
	--restrictAllelesTo BIALLELIC \
	-select "AN > 134" \
	-o ${Canids_PostFilter}

sleep 200
