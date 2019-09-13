
#Concatenate all indviduals by chromosome
#for f in {1..38}
#do 

#awk '{print FILENAME"\t"$0}' annotatedGTwithVEP_*_Chr"$f"_addROHAnnot_GERPScore.txt | grep -v "CHROM" >> AllChroms/allIndivs_Chr"$f"_annotatedGTwithVEP.txt

#done

#Add headerline back
#sed -i 1i'ID\tCHROM\tPOS\tGT\tANNOT\tIMPACT\tSIFT\twithinROH_10Kb\twithinROH_1Mb\tGERPScore' AllChroms/allIndivs_Chr*_annotatedGTwithVEP.txt

