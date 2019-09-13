

for f in {1..38}
do

#zip files
#bgzip -c masterFile_allCanidsN75_ANmt135_Chr"$f".vcf > masterFile_allCanidsN75_ANmt135_Chr"$f".vcf.gz


#index with tabix
#tabix -p vcf masterFile_allCanidsN75_ANmt135_Chr"$f".vcf.gz


#delete unzipped files
#rm -f masterFile_allCanidsN75_ANmt135_Chr"$f".vcf.idx
#rm -f masterFile_allCanidsN75_ANmt135_Chr"$f".vcf


done
