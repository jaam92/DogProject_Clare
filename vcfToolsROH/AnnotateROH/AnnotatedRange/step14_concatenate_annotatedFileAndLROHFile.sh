

#First add headerline 
#sed -i 1i'AUTO_START\tAUTO_END\tGoodSitesInROH\tROH_length\tProp_Covered' *_goodSitesWithinROH.txt

#Next grab info from original ROH file
#for i in {AW,BC,EW,IR,LB,PG,TM}; do for f in {1..38}; do paste "$i"_Chr"$f"_goodSitesWithinROH.txt ../../"$i"_Chr"$f".LROH | awk '{print $6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' | grep -v "AUTO_START" >> "$i"_allChroms_ROHQC.LROH;done;done

#Add new header line
#sed -i  1i'CHROM\tAUTO_START\tAUTO_END\tGoodSitesInROH\tROH_length\tProp_Covered\tMIN_START\tMAX_END\tN_VARIANTS_BETWEEN_MAX_BOUNDARIES\tN_MISMATCHES\tINDV' *_allChroms_ROHQC.LROH *_allChroms_ROHQC.LROH

