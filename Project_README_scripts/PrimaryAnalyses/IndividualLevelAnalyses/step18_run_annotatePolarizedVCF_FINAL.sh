#Generate Final counts of annotations
. /u/local/Modules/default/init/modules.sh
module load R/3.4.2

Rscript annotatePolarizedVCF_FINAL.R
