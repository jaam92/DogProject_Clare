#grab the CV error from job output file
grep -w "CV error" ADMIXTURE.o579487 | awk '{print $3"\t"$4}' | sed -e 's/(K=//g' | sed -e s'/)://g' | sed -e 1i'K_subpops\tCV_error' > outputFiles/CrossValidationError.out
