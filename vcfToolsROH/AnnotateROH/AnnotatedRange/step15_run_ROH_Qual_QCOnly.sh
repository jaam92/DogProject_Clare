#Generate Final ROH that pass our filters and our 10kb in length at a minimum
. /u/local/Modules/default/init/modules.sh
module load R/3.4.2

Rscript ROH_Qual_QCOnly.R
