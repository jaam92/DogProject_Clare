#Generate Final counts of annotations
. /u/local/Modules/default/init/modules.sh
module load R/3.4.2

#Generate input files
#Rscript generateG12File.R

#Run nandita's scripts
Rscript RunAnalyses/RScriptToCallandRunG12.R

