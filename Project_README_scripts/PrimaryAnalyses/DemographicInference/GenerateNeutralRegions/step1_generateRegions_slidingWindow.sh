#Want to do all of this on SCRATCH DIRECTORY

. /u/local/Modules/default/init/modules.sh
module load bcftools
#module load python/2.7

for f in {1..38} 
do 

#Step 1


###USE THESE FILES FOR SLIDING WINDOW#####
#Pull the -s will pull the 4 specified individuals; -R only grab sites within Tanya's Neutral Regions
bcftools view -s EW4,EW5,EW6,EW7 -R /u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/CanFam3.1_ExonIntronNeutral/NeutralRegions/chr"$f"_outside.ensemble_rm.phastConsEl_cutoff_0.4.bed /u/home/c/cmarsden/project-kirk-bigdata/FINAL_VCFS_May2018_INCrels/final_filtered_vcf_allsites_6files/TenEW/6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenEW_jointcalled_chr"$f".vcf.gz -o $SCRATCH/EWSamps_Chr"$f"_NeutralRegions_4slidingWindow.vcf

####USE THESE FILES FOR COMPUTING PI AND S#####
#Pull the -s will pull the 4 specified individuals; -R only grab sites within Tanya's Neutral Regions; -m2 -M2 -v snps specifies biallelic sites; --min-ac 1:minor must see the minor allele (REF or ALT) at least once
#bcftools view -s EW4,EW5,EW6,EW7 -R /u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/CanFam3.1_ExonIntronNeutral/NeutralRegions/chr"$f"_outside.ensemble_rm.phastConsEl_cutoff_0.4.bed /u/home/c/cmarsden/project-kirk-bigdata/FINAL_VCFS_May2018_INCrels/final_filtered_vcf_allsites_6files/TenEW/6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenEW_jointcalled_chr"$f".vcf.gz -m2 -M2 -v snps --min-ac 1:minor -o $SCRATCH/EWSamps_Chr"$f"_NeutralRegions.vcf

#Step2
#bgzip and tabix index

#bgzip -c EWSamps_Chr"$f"_NeutralRegions.vcf > EWSamps_Chr"$f"_NeutralRegions.vcf.gz
#bgzip -c EWSamps_Chr"$f"_NeutralRegions_4slidingWindow.vcf > EWSamps_Chr"$f"_NeutralRegions_4slidingWindow.vcf.gz
#tabix -p vcf EWSamps_Chr"$f"_NeutralRegions.vcf.gz
#tabix -p vcf EWSamps_Chr"$f"_NeutralRegions_4slidingWindow.vcf.gz

#Step 3
#Get hets per 1kb

#python /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/scripts/SlidingWindowHet_Unphased.jar.ab.jam.py --vcf EWSamps_Chr"$f"_NeutralRegions_4slidingWindow.vcf.gz --window_size 1000 --step_size 1000 --chromNum 'chr'$f''

done


