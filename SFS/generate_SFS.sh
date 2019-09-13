##Took about a half hour to run locally

source /u/local/Modules/default/init/modules.sh
module load R/3.4.2 

#generate counts
#Rscript PolarizedVCFToSFS.R
Rscript PolarizedVCFToSFS_Neutral.R

