

#mkdir $SCRATCH/EWvsTM_windowed

#load vcftools
. /u/local/Modules/default/init/modules.sh
module load vcftools

#windowed fst with 1Kb window and 1Kb step
#vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/EW_indivList.txt --weir-fst-pop /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/TM_indivList.txt --fst-window-size 1000 --fst-window-step 1000 --out EW_vs_TM

#windowed pi in TM with 1Kb window and 1Kb step
vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --keep /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/TM_indivList.txt --window-pi 1000 --window-pi-step 1000 --out TM_windowedPI_1Kbwin_1Kbstep 

#windowed pi in EW with 1Kb window and 1Kb step
vcftools --gzvcf /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/masterFile_allChroms.vcf.gz --keep /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/EW_indivList_downsample_N9.txt --window-pi 1000 --window-pi-step 1000 --out EW_N9_windowedPI_1Kbwin_1Kbstep
