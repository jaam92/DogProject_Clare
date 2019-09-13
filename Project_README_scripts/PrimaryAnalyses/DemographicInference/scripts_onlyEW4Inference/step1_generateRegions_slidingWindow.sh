
for f in {1..38} 
do 

#Step 1
#Pull EW4 and only grab sites within Tanya's Neutral Regions

#. /u/local/Modules/default/init/modules.sh
#module load bcftools
#bcftools view -s EW4 -R /u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/CanFam3.1_ExonIntronNeutral/NeutralRegions/chr"$f"_outside.ensemble_rm.phastConsEl_cutoff_0.4.bed /u/home/c/cmarsden/project-kirk-bigdata/FINAL_VCFS_May2018_INCrels/final_filtered_vcf_allsites_6files/TenEW/6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_TenEW_jointcalled_chr"$f".vcf.gz -o EW4_Chr"$f"_NeutralRegions.vcf

#Step2
#bgzip and tabix index

#bgzip -c EW4_Chr"$f"_NeutralRegions.vcf > EW4_Chr"$f"_NeutralRegions.vcf.gz
#tabix -p vcf EW4_Chr"$f"_NeutralRegions.vcf.gz

#Step 3
#Get hets per 1kb
#. /u/local/Modules/default/init/modules.sh
#module load python
#python /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/scripts/SlidingWindowHet_Unphased.jar.ab.jam.py --vcf EW4_Chr"$f"_NeutralRegions.vcf.gz --window_size 1000 --step_size 1000 --chromNum 'chr'$f''

#Step 4 
#make file for Diegos scripts
#cat EW4_Chr"$f"_NeutralRegions_het_1000win_1000step.txt | grep -v "chromo" >> /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/RunABC/InputData/EW4_allChroms_NeutralRegions_het_1000win_1000step.txt

done

