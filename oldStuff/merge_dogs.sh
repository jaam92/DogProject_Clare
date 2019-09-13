## Merge all individuals
##deleted the merged files after master files were made

CHROM=$1 ##chromosome number will be given when I qsub script

##Set working directories and input/output files 
vcfAW='/u/home/c/cmarsden/nobackup-kirk/clare/analyses/FINAL_VCFS/final_filtered_vcf_allsites_6files/TenAW/6_bespoke_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenAW_jointcalled_chr'${CHROM}'.vcf.gz'
vcfEW='/u/home/c/cmarsden/nobackup-kirk/clare/analyses/FINAL_VCFS/final_filtered_vcf_allsites_6files/TenEW/6_bespoke_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenEW_jointcalled_chr'${CHROM}'.vcf.gz'
vcfIR='/u/home/c/cmarsden/nobackup-kirk/clare/analyses/FINAL_VCFS/final_filtered_vcf_allsites_6files/TenIR/6_bespoke_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenIR_jointcalled_chr'${CHROM}'.vcf.gz'
vcfGS='/u/home/c/cmarsden/nobackup-kirk/clare/analyses/FINAL_VCFS/final_filtered_vcf_allsites_6files/TenGS/6_bespoke_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenGS_jointcalled_chr'${CHROM}'.vcf.gz'
vcfPG='/u/home/c/cmarsden/nobackup-kirk/clare/analyses/FINAL_VCFS/final_filtered_vcf_allsites_6files/TenPG/6_bespoke_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenPG_jointcalled_chr'${CHROM}'.vcf.gz'
vcfTM='/u/home/c/cmarsden/nobackup-kirk/clare/analyses/FINAL_VCFS/final_filtered_vcf_allsites_6files/TenTM/6_bespoke_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenTM_jointcalled_chr'${CHROM}'.vcf.gz'
OUTFILE='/u/flashscratch/j/jmooney3/DogsROH/mergedVCFs/6_bespokefiltered_v4_mergeGaTKfiltered_varnonvar_allSpecies_jointcalled_chr'${CHROM}'.vcf.gz'

gatk_dir='/u/home/c/cmarsden/project-klohmueldata/clares_data/software/gatk_35'
ref_genome='/u/home/c/cmarsden/project-klohmueldata/clares_data/canids/reference/canFam3/canFam3.fa'

##load java and set parameters
. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77
MEM=4G
	
### merge chromsomes from all pops
java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -T CombineVariants -R ${ref_genome} \
	--variant ${vcfAW} \
	--variant ${vcfEW} \
	--variant ${vcfIR} \
	--variant ${vcfGS} \
	 --variant ${vcfPG} \
	 --variant ${vcfTM} \
	-o ${OUTFILE} \
	--genotypemergeoption REQUIRE_UNIQUE



sleep 200
