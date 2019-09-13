
#for f in {1..38}
#do 

#grab vcf position and the reference allele then remove header
#awk '{print $3"\t"$6}' /u/home/c/cmarsden/nobackup-kirk/clare/analyses/step17_alternatereferencemaker/WD_chr"$f"_CHROM_POS_allele.txt | grep  -v "bed_end1based" > /u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/WildDog/WD_chr"$f"_CHROM_POS_allele.txt

#Grab sites in master file
#grep -v "#" /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/master_vcfClare/reformatted_masterFile_allCanidsN75_ANmt135_Chr"$f".vcf | awk '{print $2"\t"$3"\t"$4}' > /u/home/j/jmooney3/klohmueldata/jazlyn_data/Reference_Files/WildDog/CanFam3.1_ClareProj/masterFile_CanFam3.1_chr"$f"_REFALT.txt

#done

#Filter to sites in both files
#. /u/local/Modules/default/init/modules.sh
#module load R/3.4.2

#Rscript IntersectWDandMF.R


