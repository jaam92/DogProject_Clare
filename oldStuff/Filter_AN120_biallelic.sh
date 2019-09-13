

CHROM=$1 ##chromosome number will be given when I qsub script

##Set working directories and input/output files 
gatk_dir='/u/home/c/cmarsden/project-klohmueldata/clares_data/software/gatk_35'
ref_genome='/u/home/c/cmarsden/project-klohmueldata/clares_data/canids/reference/canFam3/canFam3.fa'

Mergedvcf_in='/u/flashscratch/j/jmooney3/DogsROH/mergedVCFs/6_bespokefiltered_v4_mergeGaTKfiltered_varnonvar_allSpecies_jointcalled_chr'${CHROM}'.vcf.gz'
Mastervcf_out='/u/flashscratch/j/jmooney3/DogsROH/masterVCFs/allSpecies_MasterFile_varnonvar_jointcalled_biallelic_chr'${CHROM}'.vcf.gz'

##load java and set parameters
. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77
MEM=4G

##Filter to biallelic SNPs
##Filter to sites in at least 90% of individuals (80 total so 72 is 90% covered) 
java -Xmx${MEM} -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
	-T SelectVariants \
	--variant ${Mergedvcf_in} \
	-selectType SNP \
	--restrictAllelesTo BIALLELIC \
	-select "AN > 107" \
	-o ${Mastervcf_out}

sleep 200
