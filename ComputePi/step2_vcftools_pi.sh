## qsub -cwd -V -N calcPi -l highp,time=10:00:00,h_data=1G -t 1-38:1 -M eplau -m a vcftools_pi.sh




for i in {TM,IR,LB,EW,BC,AW,PG}
do
#set chroms as num of tasks
CHROM=${SGE_TASK_ID}

#calculate relatedness using king statistic
. /u/local/Modules/default/init/modules.sh
module load vcftools

vcftools --gzvcf /u/flashscratch/j/jmooney3/AllSitesFiles/SubsetN6/"$i"'_allSitesN6_ANmt135_Chr'${CHROM}'.vcf.gz' --keep /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/IndivFiles/"$i"_indivList_downsample_N6.txt --site-pi --chr 'chr'${CHROM} --out /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/ComputePi/"$i"'_downsampledN6_chr'${CHROM}


done



