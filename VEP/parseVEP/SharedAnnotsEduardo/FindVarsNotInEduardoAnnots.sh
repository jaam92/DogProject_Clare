#Grab Vars without a snpEff annot
#awk '{print $1"_"$2, $3, $10}' ../FinalAnnot_VEPSIFTSnpEff_allChroms.txt | grep -w "NA"  | sed -e 1i'CHROM_POS\tANNOT\tSNPEff' > VarsWithoutSnpEffAnnot.txt

#Get Counts of Vars without snpEff annot
#awk '{print $3, $10}' ../FinalAnnot_VEPSIFTSnpEff_allChroms.txt | grep -w "NA" | awk '{print $1}' | sort | uniq -c > CountsVarsWithoutSnpEffAnnot.txt

#Pull Eduardos annotations for those sites with only 3 transcripts
#awk '$2 == 3 {print $0}'  /u/home/j/jmooney3/klohmueldata/cegamori/0_and_4_deg_sites/coord_annot.list > 3annotsEduardo.txt

#Compare the two 
#. /u/local/Modules/default/init/modules.sh
#module load R/3.4.2
#Rscript CompareAnnots.R 

