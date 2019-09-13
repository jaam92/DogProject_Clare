#Reformat annots for Clare's script
#LOF,SY,NS is what clares script expects
#sed -e 's/missense_variant,splice_region_variant/missense_variant/g' FinalAnnot_VEPSIFTSnpEff_allChroms.txt | sed -e 's/splice_donor_variant,missense_variant/missense_variant/g' | sed -e 's/splice_region_variant,synonymous_variant/synonymous_variant/g' | sed -e 's/splice_region_variant,synonymous_variant/synonymous_variant/g' | sed -e 's/stop_retained_variant/synonymous_variant/g' | sed -e 's/stop_gained,splice_region_variant/stop_gained/g' | sed -e 's/stop_gained,start_lost/stop_gained/g' | sed -e 's/stop_lost/stop_gained/g' | sed -e 's/start_lost/stop_gained/g' | sed -e 's/start_lost,synonymous_variant/stop_gained/g' | sed -e 's/stop_gained,start_lost/stop_gained/g' | sed -e 's/missense_variant/NS/g' | sed -e 's/synonymous_variant/SY/g' | sed -e 's/stop_gained/LOF/g' | sed -e 's/LOF,SY/LOF/g' | grep -v "splice" > FinalAnnot_VEPSIFTSnpEff_allChroms_updatedPTVSplice.txt

#Change one of the annotations with a weird gene name gives error for column count downstream
#sed -i 's/C22H13orf46 (NCBI)/C22H13orf46/g' FinalAnnot_VEPSIFTSnpEff_allChroms_updatedPTVSplice.txt


#Split annotations by chromosome

#for f in {1..38}
#do

#grep -w "chr$f" FinalAnnot_VEPSIFTSnpEff_allChroms_updatedPTVSplice.txt > parsedVEP_chr"$f"_codingAnnots.txt

#done 
