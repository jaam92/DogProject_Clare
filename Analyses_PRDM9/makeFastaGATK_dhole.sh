#Dhole01 comes from VCF from this paper:Whole genome sequencing of canids reveals genomic regions under selection and variants influencing morphology (Plassais et al. Nat Comms 2019)

#directories
gatk_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/GATK/gatk_35'
ref_genome='/u/home/c/cmarsden/project-klohmueldata/clares_data/canids/reference/canFam3/canFam3.fa'

#Load java
. /u/local/Modules/default/init/modules.sh
module load java


#Step 1 pull data from PRDM9 region and pull a single individual
java -Xmx4G -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
        -T SelectVariants \
        --variant /u/flashscratch/j/jmooney3/SRZ189891_722g.990.SNP.INDEL.chrAll.vcf.gz \
        -sn Dhole01 \
        -L chr5:63559618-63626468 \
        -o Dhole01_Plassais2019_NatComm_PRDM9.vcf

#Step 2 make our fasta file for that individual
java -Xmx4G -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
        -T FastaAlternateReferenceMaker \
        --variant Dhole01_Plassais2019_NatComm_PRDM9.vcf \
        -L chr5:63559618-63626468 \
        -o Dhole01_Plassais2019_NatComm_PRDM9.fa

#redo header
#sed -i 's/>1 chr5:63559618/>dhole/g' Dhole01_Plassais2019_NatComm_PRDM9.fa

#Didn't find any sites in GAPDH region for dhole 
#Step 1 pull data from GAPDH region and pull a single individual
#java -Xmx4G -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
#        -T SelectVariants \
#        --variant /u/flashscratch/j/jmooney3/SRZ189891_722g.990.SNP.INDEL.chrAll.vcf.gz \
#        -sn Dhole01 \
#        -L chr27:38468698-38472776 \
#        -o Dhole01_Plassais2019_NatComm_GAPDH.vcf

#Step 2 make our fasta file for that individual
#java -Xmx4G -jar ${gatk_dir}/GenomeAnalysisTK.jar -R ${ref_genome} \
#        -T FastaAlternateReferenceMaker \
#        --variant Dhole01_Plassais2019_NatComm_GAPDH.vcf \
#        -L chr27:38468698-38472776 \
#        -o Dhole01_Plassais2019_NatComm_GAPDH.fa


