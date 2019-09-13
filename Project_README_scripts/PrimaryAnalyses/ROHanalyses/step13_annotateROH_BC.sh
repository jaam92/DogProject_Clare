##qsub -cwd -V -N ROHQual -l highp,time=10:00:00,h_data=5G -t 1-38:1 -M eplau -m ea annotateROH_BC.sh

#Load python
. /u/local/Modules/default/init/modules.sh
module load python/2.7

#Set variables
CHROM=${SGE_TASK_ID}
in_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/vcfToolsROH/AnnotateROH'
out_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/vcfToolsROH/AnnotateROH/AnnotatedRange'
site_dir='/u/flashscratch/j/jmooney3/allSitesPerSpecies'

#Run script

        ROHrange_in=${in_dir}/'BC_Chr'$CHROM'_ROHRanges.txt'

        sites_in=${site_dir}/'TenBC_chr'$CHROM'_allSitesPosOnly.txt'

        count_out=${out_dir}/'BC_Chr'$CHROM'_goodSitesWithinROH.txt'

        python ROHQualCounts_v2.py ${ROHrange_in} ${sites_in} 1 ${count_out}



sleep 180
