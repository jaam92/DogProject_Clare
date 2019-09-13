for f in {1..38} 
do 

#remove stuff above header line
grep -v "##" annotatedVEPandSIFT_masterFile_allCanidsN75_ANmt135_Chr"$f".txt  > parseVEP/masterFile_allCanidsN75_ANmt135_Chr"$f"_fixHeader.txt


done
