

##qsub -cwd -V -N slidingWindows -l highp,time=00:20:00,h_data=10G -t 1-38:1 -M eplau -m a runSlidingWindow.sh
##CHROM=${SGE_TASK_ID}

for CHROM in {1..38}; 
do 

inputDir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare'
master_VCF=${inputDir}/'masterFile_allCanidsN75_ANmt135_Chr'${CHROM}'.vcf.gz'


#load python
. /u/local/Modules/default/init/modules.sh
module load python

#command to run script
#python SlidingWindowHet_Unphased.jar.ab.jam.py --vcf ${master_VCF} --window_size 1000000 --step_size 1000000 --chromNum 'chr'${CHROM}
python SlidingWindowHet_Unphased.jar.ab.jam.py --vcf ${master_VCF} --window_size 10000 --step_size 10000 --chromNum 'chr'${CHROM}
done


