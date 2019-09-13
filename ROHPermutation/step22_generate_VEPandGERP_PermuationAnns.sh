#Step 1: Generate annotation files for permutation
#. /u/local/Modules/default/init/modules.sh
#module load R/3.4.2
#Rscript AnnotationFiles/MakePermSnpFile.R

#Step 2: Split annotation file by chromosome
#for f in {1..38} 
#do 

#Split VEP
#grep -w "chr$f" AnnotationFiles/allChromsVEPAnnotation.txt | sed -e 's/chr//g' > chr"$f"_FinalVEPAnnots.txt; 

#SplitGERP
#grep -w "chr$f" AnnotationFiles/allChromsGERPScoreAnnotation.txt | awk '{print $1"\t"$2"\t"$5}' | sed -e 's/chr//g' > chr"$f"_FinalGERPAnnots.txt

#done
