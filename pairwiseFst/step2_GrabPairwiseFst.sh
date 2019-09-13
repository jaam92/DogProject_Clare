
#grab compairson groups
#grep "\\-out" /u/flashscratch/j/jmooney3/pairwiseFst/pwFst.out > comparison

#grab fst
#grep "weighted" /u/flashscratch/j/jmooney3/pairwiseFst/pwFst.out > Weighted_Fst

#Make comparison file
#paste comparison Weighted_Fst | awk '{print $2"\t"$9}' > pairwiseFst.out
