#for i in {AW10,AW11,AW12,AW13,AW15,AW16,AW17,AW18,AW19,AW20,AW21,AW22,AW23,AW24,BC1,BC3,BC4,BC5,BC6,BC7,BC9,EW1,EW10,EW2,EW3,EW4,EW5,EW6,EW7,EW8,EW9,IR1,IR10,IR2,IR3,IR4,IR5,IR6,IR7,IR8,IR9,LB1,LB11,LB2,LB3,LB4,LB5,LB6,LB7,LB8,LB9,PG1,PG10,PG11,PG12,PG13,PG14,PG16,PG2,PG3,PG4,PG5,PG6,PG7,PG8,PG9,TM1,TM10,TM2,TM3,TM5,TM6,TM7,TM8,TM9}
#do

#for f in {1..38}
#do

#cat annotatedGTwithVEP_"$i"_Chr"$f"_addROHAnnot_GERPScore.txt | grep -v "CHROM" >> AllChroms/annotatedGTwithVEP_"$i"_allChroms.txt

#done
#done

#Add headerline back
#sed -i 1i'CHROM\tPOS\tGT\tANNOT\tIMPACT\tSIFT\twithinROH_10Kb\twithinROH_1Mb\tGERPScore' AllChroms/annotatedGTwithVEP_*_allChroms.txt
