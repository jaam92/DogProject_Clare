#Step 4 
#make file for Diegos scripts

#for f in {1..38}
#do
#cat EWSamps_Chr"$f"_NeutralRegions_het_1000win_1000step.txt | grep -v "chromo" >> EWSamps_allChroms_NeutralRegions_het_1000win_1000step_noHeader.txt
#done

#Step 5
#put header back

#head -n 1 EWSamps_Chr38_NeutralRegions_het_1000win_1000step.txt > header.txt
#cat header.txt EWSamps_allChroms_NeutralRegions_het_1000win_1000step_noHeader.txt > EWSamps_allChroms_NeutralRegions_het_1000win_1000step_Header.txt

#Step 6 
#compute summary statistics files

#. /u/local/Modules/default/init/modules.sh
#module load vcftools

#for f in {1..38}
#do
#Compute pi per site
#vcftools --gzvcf $SCRATCH/EWSamps_Chr"$f"_NeutralRegions.vcf.gz --chr chr"$f" --out $SCRATCH/chr"$f" --site-pi 

#Compute S per 10 bp window
#vcftools --gzvcf $SCRATCH/EWSamps_Chr"$f"_NeutralRegions.vcf.gz --chr chr"$f" --out $SCRATCH/chr"$f"_SegSites --SNPdensity 10
#done

#Step 7
#Concatenate summary statistics across chromosomes remove sites/windows with 0 values

#for f in {1..38}
#do
#grep -v "CHROM" $SCRATCH/chr"$f".sites.pi | awk '$3 != 0 {print $0}' >> allChroms.sites.pi
#grep -v "CHROM" $SCRATCH/chr"$f"_SegSites.snpden | awk '$3 != 0 {print $1"\t"$2"\t"$3}' >> allChroms.snpden 
#done

#Step 8 
#add header to file

#sed -i 1i'CHROM\tPOS\tPI' allChroms.sites.pi 
#sed -i 1i'CHROM\tPOS\tSNPCOUNT' allChroms.snpden


#Step 9
#make input file for ABC

. /u/local/Modules/default/init/modules.sh
module load R/3.4.2
Rscript generateABCInput_4samps.R


